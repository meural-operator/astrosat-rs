//! Unit tests covering every stage of the AstroSat LAXPC Level-1 -> Level-2 pipeline.

use astrosat_core::gti::{GtiInterval, GtiSet};
use astrosat_core::mkf::MkfTable;
use astrosat_core::tct::{LaxpcTct, TctTable};
use astrosat_fits::{validate_ogip, FitsFile, OgipMeta};
use astrosat_laxpc::*;
use std::collections::HashMap;
use std::path::Path;
use tempfile::tempdir;

/// Helper to construct a valid synthetic 2242-byte EA frame with configurable header & payload words.
fn make_ea_frame(
    l1_time: f64,
    unit_id: u8,
    frame_counter: u16,
    stbg_40: u64,
    ie1: i16,
    words: &[[u8; 5]],
    add_eof: bool,
) -> Vec<u8> {
    let mut raw = vec![0u8; L1_FRAME_SIZE];
    raw[0..8].copy_from_slice(&l1_time.to_be_bytes());
    raw[44] = FRAME_SYNC_BYTE;
    raw[45] = unit_id;
    raw[46] = MODE_EA_BYTE;
    let s40_be = stbg_40.to_be_bytes();
    raw[49..54].copy_from_slice(&s40_be[3..8]);
    raw[57] = 0; // valid status flag
    raw[58..60].copy_from_slice(&frame_counter.to_be_bytes());

    let mut pos = 60;
    for w in words {
        raw[pos..pos + 5].copy_from_slice(w);
        pos += 5;
    }
    if add_eof && pos < 2092 {
        raw[pos] = END_OF_FRAME_MARKER;
    }
    raw[2092..2094].copy_from_slice(&ie1.to_be_bytes());
    raw
}

/// Helper to encode a 5-byte EA event word matching `RawEaWord::unpack`.
fn encode_ea_word(fine_ticks: u8, ia1: u8, ia2: u8, ic1: u16, ic2: u16, k1: u8, k2: u8) -> [u8; 5] {
    let b0 = (ia1 & 0x0F) | ((ia2 & 0x0F) << 4);
    let b1 = fine_ticks.wrapping_add(2);
    let b2 = ((ic1 >> 2) & 0xFF) as u8;
    let b3 = ((ic2 >> 2) & 0xFF) as u8;
    let l1 = (ic1 & 0x03) as u8;
    let l2 = (ic2 & 0x03) as u8;
    let b4 = (k1 & 0x01) | (l1 << 2) | ((k2 & 0x01) << 4) | (l2 << 6);
    [b0, b1, b2, b3, b4]
}

#[test]
fn stage1_frame_sorting_rollover_dedup_stray_and_gaps() {
    // Out-of-order frames crossing 16-bit wrap + duplicate + stray + corrupted header + gap
    // Valid frames: 100, 102 (out of order), 101, 102 (duplicate), 20000 (stray!), 103, 105 (gap at 104)
    let mut frames_raw = vec![
        make_ea_frame(1000.0, 1, 100, 1_000_000, 0, &[], true),
        make_ea_frame(1002.0, 1, 102, 1_020_000, 0, &[], true),
        make_ea_frame(1001.0, 1, 101, 1_010_000, 0, &[], true),
        make_ea_frame(1002.01, 1, 102, 1_020_000, 0, &[], true), // duplicate replacement
        make_ea_frame(1002.5, 1, 20_000, 1_025_000, 0, &[], true), // stray jump
        make_ea_frame(1003.0, 1, 103, 1_030_000, 0, &[], true),
        make_ea_frame(1005.0, 1, 105, 1_050_000, 0, &[], true), // missing 104
    ];

    // Corrupted header (wrong sync)
    let mut bad_sync = make_ea_frame(1006.0, 1, 106, 1_060_000, 0, &[], true);
    bad_sync[44] = 0xAA;
    frames_raw.push(bad_sync);

    // Corrupted RS error flag (ie1 = 2 > 1)
    frames_raw.push(make_ea_frame(1006.0, 1, 106, 1_060_000, 2, &[], true));

    // Wrong unit (unit 2 when requesting unit 1)
    frames_raw.push(make_ea_frame(1006.0, 2, 106, 1_060_000, 0, &[], true));

    let slices: Vec<&[u8]> = frames_raw.iter().map(|v| v.as_slice()).collect();
    let (ordered, diag) =
        sort_and_dedup_frames(slices, Some(LaxpcUnit::Lxp10), Some(FrameMode::Ea)).unwrap();

    let counters: Vec<u16> = ordered.iter().map(|f| f.header.frame_counter).collect();
    assert_eq!(counters, vec![100, 101, 102, 103, 105]);
    assert_eq!(diag.rejected_header, 3);
    assert_eq!(diag.rejected_stray, 1);
    assert_eq!(diag.duplicates_replaced, 1);
    assert_eq!(diag.single_frame_gaps, 1);
    assert_eq!(diag.total_missing_frames, 1);
    assert_eq!(diag.gaps.len(), 1);
    assert_eq!(diag.gaps[0].prev_counter, 103);
    assert_eq!(diag.gaps[0].next_counter, 105);
}

#[test]
fn stage2_ea_5byte_unpacking_and_double_event_logic() {
    // 1. Single event (ic2 = 0)
    let w_single = encode_ea_word(42, 3, 0, 250, 0, 0, 0);
    let raw = RawEaWord::unpack(w_single);
    assert_eq!(raw.fine_ticks, 42);
    assert_eq!(raw.ia1, 3);
    assert_eq!(raw.ic1, 250);
    assert_eq!(raw.ic2, 0);
    assert_eq!(raw.resolve(), Some((250, 3, EventKind::Single)));

    // 2. Xe K-escape double event: k1=1, ic1 <= ic2 -> assigns ia2
    let w_dbl1 = encode_ea_word(10, 1, 4, 240, 310, 1, 0);
    assert_eq!(
        RawEaWord::unpack(w_dbl1).resolve(),
        Some((550, 4, EventKind::Double))
    );

    // 3. Xe K-escape double event: k1=1, k2=0, ic1 > ic2 > 1 -> assigns ia2 (`laxpcl1.f:2312-2313`)
    let w_dbl2 = encode_ea_word(11, 2, 5, 300, 100, 1, 0);
    assert_eq!(
        RawEaWord::unpack(w_dbl2).resolve(),
        Some((400, 5, EventKind::Double))
    );

    // 4. Xe K-escape double event: k1=0, k2=1 -> assigns ia1
    let w_dbl3 = encode_ea_word(12, 2, 5, 400, 240, 0, 1);
    assert_eq!(
        RawEaWord::unpack(w_dbl3).resolve(),
        Some((640, 2, EventKind::Double))
    );

    // 5. Double event PHA saturation clamped at 1023
    let w_sat = encode_ea_word(13, 1, 2, 700, 500, 0, 1);
    assert_eq!(
        RawEaWord::unpack(w_sat).resolve(),
        Some((1023, 1, EventKind::Double))
    );

    // 6. Non-K double event (k1=0, k2=0, ic2>0) -> rejected
    let w_non_k = encode_ea_word(14, 1, 2, 200, 150, 0, 0);
    assert_eq!(RawEaWord::unpack(w_non_k).resolve(), None);

    // 7. Veto double event (ia2 = 8 >= 8) -> rejected even with k1=1
    let w_veto = encode_ea_word(15, 1, 8, 240, 200, 1, 0);
    assert_eq!(RawEaWord::unpack(w_veto).resolve(), None);
}

#[test]
fn stage2_ea_frame_decode_coarse_markers_and_gaps() {
    // Build 2 consecutive frames with a 0.25 s gap between them
    let stbg1: u64 = 10_000_000; // 100.0 s in 100 kHz ticks
    let l1_t1 = 1.8121823219164e8 + stbg1 as f64 / STBG_CLOCK_HZ;

    // Coarse marker at 10_000_000 >> 8 = 39062 -> m4=0, m3=152, m2=150
    let c1 = (stbg1 >> 8) as u32;
    let ef1 = [
        COARSE_TIME_MARKER,
        ((c1 >> 16) & 0xFF) as u8,
        ((c1 >> 8) & 0xFF) as u8,
        (c1 & 0xFF) as u8,
        0,
    ];
    let w1 = encode_ea_word(10, 1, 0, 120, 0, 0, 0);
    let w2 = encode_ea_word(50, 2, 0, 220, 0, 0, 0);
    let w3 = encode_ea_word(90, 3, 0, 320, 0, 0, 0);

    let f1_raw = make_ea_frame(l1_t1, 1, 1, stbg1, 0, &[ef1, w1, w2, w3], false);

    // Second frame 0.25 s later (+25,000 ticks)
    let stbg2: u64 = stbg1 + 25_000;
    let l1_t2 = 1.8121823219164e8 + stbg2 as f64 / STBG_CLOCK_HZ;
    let c2 = (stbg2 >> 8) as u32;
    let ef2 = [
        COARSE_TIME_MARKER,
        ((c2 >> 16) & 0xFF) as u8,
        ((c2 >> 8) & 0xFF) as u8,
        (c2 & 0xFF) as u8,
        0,
    ];
    let w4 = encode_ea_word(20, 1, 0, 150, 0, 0, 0);
    let f2_raw = make_ea_frame(l1_t2, 1, 2, stbg2, 0, &[ef2, w4], false);

    let (ordered, _) = sort_and_dedup_frames(
        [f1_raw.as_slice(), f2_raw.as_slice()],
        Some(LaxpcUnit::Lxp10),
        Some(FrameMode::Ea),
    )
    .unwrap();

    let stream = decode_ea_orbits(&[ordered], LaxpcUnit::Lxp10);
    assert_eq!(stream.events.len(), 4);
    assert_eq!(stream.events[0].pha, 120);
    assert_eq!(stream.events[1].pha, 220);
    assert_eq!(stream.events[2].pha, 320);
    assert_eq!(stream.events[3].pha, 150);
    // Internal frame gap >= 50 ms should be detected
    assert_eq!(stream.frame_gaps.len(), 1);
    let gap_dur = stream.frame_gaps[0].1 - stream.frame_gaps[0].0;
    assert!((gap_dur - 0.25).abs() < 0.01, "expected ~0.25 s gap, got {gap_dur}");
}

#[test]
fn stage3_bbc_subframe_decoding_and_quartiles() {
    let mut raw = vec![0u8; L1_FRAME_SIZE];
    raw[44] = FRAME_SYNC_BYTE;
    raw[45] = 1;
    raw[46] = MODE_BB_BYTE;

    // Populate 31 subframes with 128 ms cadence starting at t_ins = 200_000.0 s
    let t0_ticks: u64 = 20_000_000_000;
    for j in 1..=31usize {
        let base = j * 64;
        let ticks = t0_ticks + ((j - 1) as u64) * 12_800; // 0.128 s per subframe
        let tb = ticks.to_be_bytes();
        raw[base + 47..base + 54].copy_from_slice(&tb[1..8]);
        // LLD=500, JAN=300, LDU=110 (220 per pair -> ~859 cps -> quartile 2), IGE=40
        raw[base + 56..base + 58].copy_from_slice(&500i16.to_be_bytes());
        raw[base + 58..base + 60].copy_from_slice(&300i16.to_be_bytes());
        raw[base + 60..base + 62].copy_from_slice(&110i16.to_be_bytes());
        raw[base + 64..base + 66].copy_from_slice(&40i16.to_be_bytes());
    }

    let recs = decode_bb_frame(&raw).unwrap();
    assert_eq!(recs.len(), 16);
    // First odd subframe (j=1): single subframe counts
    assert_eq!(recs[0].lld, 500);
    assert_eq!(recs[0].jan, 300);
    assert_eq!(recs[0].ldu, 110);
    assert_eq!(recs[0].ige, 40);
    assert_eq!(recs[0].total_rate_counts, 300 + 110 + 40);
    // Subsequent odd subframes (j=3, 5, ...): paired sums
    assert_eq!(recs[1].lld, 1000);
    assert_eq!(recs[1].jan, 600);
    assert_eq!(recs[1].ldu, 220);
    assert_eq!(recs[1].ige, 80);

    // Check quartile classification for LAXPC10 (bounds 750..1050, dul=75)
    assert_eq!(BbcSeries::uld_quartile_index(LaxpcUnit::Lxp10, 700.0), 0);
    assert_eq!(BbcSeries::uld_quartile_index(LaxpcUnit::Lxp10, 780.0), 1);
    assert_eq!(BbcSeries::uld_quartile_index(LaxpcUnit::Lxp10, 860.0), 2);
    assert_eq!(BbcSeries::uld_quartile_index(LaxpcUnit::Lxp10, 950.0), 3);
    assert_eq!(BbcSeries::uld_quartile_index(LaxpcUnit::Lxp10, 1020.0), 4);
    assert_eq!(BbcSeries::uld_quartile_index(LaxpcUnit::Lxp10, 1200.0), 5);
}

#[test]
fn stage4_tct_roundtrip_and_standard_gti_elv_saa_merge_drop() {
    // Construct synthetic TCT with 1 glitch row (ratio 1.05 > 1.001) and normal ratio 1.0002
    let ins = vec![1000.0, 1016.0, 1032.0, 1048.0, 1064.0, 1200.0];
    let sps = vec![
        189_900_000.0,
        189_900_016.003_2,
        189_900_033.5, // glitch! dt_sps/dt_ins = 1.093 -> rejected
        189_900_048.009_6,
        189_900_064.012_8,
        189_900_200.040_0,
    ];
    let tct = LaxpcTct::new(TctTable::from_pairs(sps, ins).unwrap()).unwrap();
    assert_eq!(tct.rejected_glitch_count(), 1);

    let t_utc = tct.interpolate(1080.5);
    let t_ins_back = utc_to_ins(&tct, t_utc);
    assert!(
        (t_ins_back - 1080.5).abs() < 1e-6,
        "utc_to_ins roundtrip error: got {t_ins_back}"
    );

    // Build 100 s of 1-s bins:
    // 0..30s: good (ELV=30, LON=100, rate=500)
    // 30..45s: Earth occulted (ELV=-5)
    // 45..48s: 3s short good island (< 5 s -> dropped!)
    // 48..65s: SAA (LON=270 > 249)
    // 65..100s: good (ELV=25, LON=80, rate=500)
    let n = 100;
    let t0_ins = 1000.0;
    let mut mkf_t = Vec::with_capacity(n);
    let mut elv = Vec::with_capacity(n);
    let mut lon = Vec::with_capacity(n);
    let mut lat = Vec::with_capacity(n);
    let coall = vec![500.0; n];

    for i in 0..n {
        let tc = t0_ins + i as f64 + 0.5;
        mkf_t.push(tct.interpolate(tc));
        lat.push(0.0);
        if (30..45).contains(&i) {
            elv.push(-5.0);
            lon.push(100.0);
        } else if (48..65).contains(&i) {
            elv.push(25.0);
            lon.push(270.0);
        } else {
            elv.push(30.0);
            lon.push(100.0);
        }
    }
    // Add a 0.2s internal frame gap inside first good interval
    let frame_gaps = vec![(1010.0, 1010.2)];

    let mut cols = HashMap::new();
    cols.insert("ELV".into(), elv);
    cols.insert("EARTHLON".into(), lon);
    cols.insert("EARTHLAT".into(), lat);
    let mkf = MkfTable::from_columns(mkf_t, cols).unwrap();

    let raw = RawBinnedCounts {
        t0_ins,
        tbin: 1.0,
        selected_counts: coall.clone(),
        all_main_counts: coall,
    };

    let cfg = GtiConfig {
        enable_diurnal_saa_extension: false,
        ..Default::default()
    };
    let res = generate_standard_gti(
        LaxpcUnit::Lxp10,
        &raw,
        &tct,
        &mkf,
        None,
        &frame_gaps,
        &cfg,
        None,
    )
    .unwrap();

    // Expect exactly 2 intervals: [1000..1030] (30s) and [1065..1100] (35s); 3s island dropped!
    assert_eq!(res.gti_ins.len(), 2, "got intervals: {:?}", res.gti_ins.intervals());
    assert!((res.gross_exposure_sec - 65.0).abs() < 1.5);
    assert!((res.frame_gap_loss_sec - 0.2).abs() < 1e-6);
    assert!((res.net_exposure_sec - (res.gross_exposure_sec - 0.2)).abs() < 1e-6);

    // User GTI intersection restricts to [t(1005), t(1020)] -> 15 s
    let user_gti = GtiSet::new(vec![GtiInterval::new(
        tct.interpolate(1005.0),
        tct.interpolate(1020.0),
    )]);
    let res_user = generate_standard_gti(
        LaxpcUnit::Lxp10,
        &raw,
        &tct,
        &mkf,
        None,
        &frame_gaps,
        &cfg,
        Some(&user_gti),
    )
    .unwrap();
    assert_eq!(res_user.gti_ins.len(), 1);
    assert!((res_user.gross_exposure_sec - 15.0).abs() < 1e-4);
}

#[test]
fn stage5_deadtime_math_smoothing_and_ogip_fits_roundtrip() {
    // Non-paralyzable deadtime at 43 us and 42 us
    let r43 = deadtime_correct_rate(2000.0, 2000.0, 43.0e-6);
    assert!((r43 - 2000.0 / (1.0 - 0.086)).abs() < 1e-9);
    let r42 = deadtime_correct_rate(2000.0, 2000.0, 42.0e-6);
    assert!((r42 - 2000.0 / (1.0 - 0.084)).abs() < 1e-9);

    // Singularity guard at >= 23000 cps returns uncorrected rate (no Inf/NaN/negative)
    let r_guard = deadtime_correct_rate(25_000.0, 25_000.0, 43.0e-6);
    assert_eq!(r_guard, 25_000.0);

    // Running 1-s mean smoothing when tbin = 0.05 s (< 0.2 s)
    let mut coall = vec![10.0; 100];
    coall[50] = 220.0; // single spike
    let smoothed = smooth_all_counts_for_deadtime(&coall, 0.05);
    assert!(smoothed[50] < 30.0 && smoothed[50] > 15.0);
    // No smoothing when tbin >= 0.2 s
    let unsmoothed = smooth_all_counts_for_deadtime(&coall, 0.5);
    assert_eq!(unsmoothed[50], 220.0);

    // OGIP Format B & Format A FITS roundtrip + validator check
    let dir = tempdir().unwrap();
    let meta = OgipMeta {
        instrume: "LAXPC".into(),
        detnam: Some("LX10".into()),
        object: Some("CRAB".into()),
        obs_id: Some("9000000252".into()),
        ra_obj: Some(83.6331),
        dec_obj: Some(22.0145),
        tstart: 189_907_200.0,
        tstop: 189_907_300.0,
        exposure: 100.0,
        ..Default::default()
    };

    let evs = vec![
        CalibratedEvent {
            time_utc: 189_907_201.5,
            channel: 45,
            anode_id: 1,
            layer: 1,
            detector_id: 1,
        },
        CalibratedEvent {
            time_utc: 189_907_202.25,
            channel: 120,
            anode_id: 3,
            layer: 2,
            detector_id: 2,
        },
    ];
    let gti_utc = vec![(189_907_200.0, 189_907_300.0)];

    let p_evn = dir.path().join("lxp1level2.evn");
    write_format_b_event_fits(
        &p_evn,
        &evs,
        &gti_utc,
        LaxpcUnit::Lxp10,
        &EventFilter::default(),
        512,
        &meta,
    )
    .unwrap();
    let f_evn = FitsFile::open(&p_evn).unwrap();
    assert!(validate_ogip(&f_evn).is_empty());
    assert_eq!(f_evn.hdu_by_name("LAXPC_MODE_EA").unwrap().as_bintable().unwrap().nrows(), 2);

    let p_fa = dir.path().join("level2.event.fits");
    write_format_a_event_fits(&p_fa, &evs, &gti_utc, &meta).unwrap();
    let f_fa = FitsFile::open(&p_fa).unwrap();
    assert!(validate_ogip(&f_fa).is_empty());
    let fa_hdu = f_fa.hdu_by_name("LAXPC_EVENTS").unwrap();
    assert_eq!(fa_hdu.header().get_i64("MJDREFI"), Some(55197));
    assert_eq!(fa_hdu.header().get_str("TIMESYS"), Some("UTC"));
    let fa_bt = fa_hdu.as_bintable().unwrap();
    assert_eq!(fa_bt.read_col_u8("LAXPC_No").unwrap(), vec![1, 2]);
    assert_eq!(fa_bt.read_col_u8("Layer").unwrap(), vec![1, 2]);

    let p_spec = dir.path().join("lxp2level2.spec");
    let rates = vec![1.5; 256];
    let errs = vec![0.1; 256];
    write_spectrum_fits(
        &p_spec,
        LaxpcUnit::Lxp20,
        &rates,
        &errs,
        100.0,
        43e-6,
        &EventFilter::default(),
        "lx20v1.0.rmf",
        "NONE",
        &meta,
    )
    .unwrap();
    let f_spec = FitsFile::open(&p_spec).unwrap();
    assert!(validate_ogip(&f_spec).is_empty());
}

#[test]
fn stage6_gain_alignment_hv_jumps_rmf_and_2d_background_model() {
    let ref_dir = std::env::var("ASTROSAT_TEST_DATA")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|_| std::path::PathBuf::from("/mnt/d/ISRO/data"))
        .join("ref_src/laxpc");
    if !ref_dir.exists() {
        eprintln!(
            "Skipping stage6 gain/background reference test (ref_src not present: {})",
            ref_dir.display()
        );
        return;
    }

    let gt = GainTable::from_file(&ref_dir.join("gaina8.dat")).unwrap();
    let hv = HvJumpTable::from_file(&ref_dir.join("lxhv.dat")).unwrap();
    assert!(gt.entries.len() > 1000);
    assert!(!hv.jumps_by_unit[0].is_empty());

    // Crab ObsID 9000000252 was observed on 2016-01-07 -> mission day ~101
    let sol10 = compute_gain_shift(
        LaxpcUnit::Lxp10,
        101.0,
        100.0,
        Some(237.3),
        &gt,
        Some(&hv),
        false,
        false,
    )
    .unwrap();
    assert!((sol10.source_peak_ch - 238.2).abs() < 3.0);
    assert!(sol10.crossed_hv_jumps.is_empty());
    assert!(sol10.recommended_rmf.starts_with("lx10"));
    assert!(sol10.recommended_rmf.ends_with(".rmf"));

    // Crossing an HV jump on LAXPC10 (day 96) between day 90 and day 101 must be flagged
    let sol_jump = compute_gain_shift(
        LaxpcUnit::Lxp10,
        101.0,
        90.0,
        None,
        &gt,
        Some(&hv),
        false,
        false,
    )
    .unwrap();
    assert!(sol_jump.crossed_hv_jumps.contains(&96.0));

    // ChannelMap energy-to-channel selection
    let cmap10 = ChannelMap::from_file(&ref_dir.join("lxp10.channel")).unwrap();
    assert_eq!(cmap10.len(), 512);
    let (p_lo, p_hi) = cmap10.energy_to_pha_1024(3.0, 80.0);
    assert!(p_lo < p_hi && p_lo > 10 && p_hi <= 1023);

    // Background spectrum gain shifting preserves total flux within ~2% for small shift
    let bkg_pha = BackgroundPha::from_file(&ref_dir.join("backlxp10jan16.pha")).unwrap();
    assert_eq!(bkg_pha.rates.len(), 512);
    let shifted = shift_background_spectrum(
        LaxpcUnit::Lxp10,
        1.5,
        1.0,
        101.0,
        &cmap10,
        &bkg_pha.rates,
    )
    .unwrap();
    let sum_orig: f64 = bkg_pha.rates.iter().map(|r| r[10]).sum();
    let sum_sh: f64 = shifted.iter().map(|r| r[10]).sum();
    assert!(
        ((sum_sh - sum_orig) / sum_orig).abs() < 0.03,
        "flux drift too large: orig={sum_orig:.2}, shifted={sum_sh:.2}"
    );

    // Real 2D B-spline + 5-harmonic diurnal background model (`backfit.lxp10mar16`)
    let model = BackgroundModel2D::from_file(&ref_dir.join("backfit.lxp10mar16")).unwrap();
    assert_eq!(model.lon_knots.len(), 40);
    assert_eq!(model.lat_knots.len(), 10);
    let r1 = model.evaluate(100.0, 0.0, 189_910_000.0);
    let r2 = model.evaluate(100.0, 0.0, 189_910_000.0 + 21_124.0); // quarter diurnal period later
    assert!(r1 > 150.0 && r1 < 350.0, "unexpected 2D bkg rate: {r1}");
    assert!((r1 - r2).abs() > 0.01, "expected diurnal modulation between {r1} and {r2}");

    // RMF recommendations across units and epochs
    assert_eq!(recommend_rmf(LaxpcUnit::Lxp10, 101.0, 239.1, false, false), "lx10v1.1.rmf");
    assert_eq!(recommend_rmf(LaxpcUnit::Lxp20, 101.0, 224.1, true, true), "lx20L1SEv1.0.rmf");
    assert!(recommend_rmf(LaxpcUnit::Lxp30, 400.0, 240.0, false, false).contains("rho"));
}

#[test]
fn audit_regressions_safety_math_and_edge_cases() {
    // 1. Corrupt negative Fortran backfit dimensions must return Err (never panic/OOM)
    let mut corrupt = Vec::new();
    let rec0: [u8; 32] = {
        let mut b = [0u8; 32];
        b[0..4].copy_from_slice(&(-5i32).to_le_bytes());
        b[4..8].copy_from_slice(&10i32.to_le_bytes());
        b[8..12].copy_from_slice(&42i32.to_le_bytes());
        b[12..16].copy_from_slice(&514i32.to_le_bytes());
        b
    };
    corrupt.extend_from_slice(&32u32.to_le_bytes());
    corrupt.extend_from_slice(&rec0);
    corrupt.extend_from_slice(&32u32.to_le_bytes());
    for _ in 0..3 {
        corrupt.extend_from_slice(&8u32.to_le_bytes());
        corrupt.extend_from_slice(&0.0f64.to_le_bytes());
        corrupt.extend_from_slice(&8u32.to_le_bytes());
    }
    assert!(BackgroundModel2D::from_fortran_bytes(&corrupt).is_err());

    // 2. Degenerate B-spline inputs must not panic
    let (l0, b0) = eval_bspline_basis(&[1.0], 4, 1.0);
    assert_eq!(l0, 0);
    assert_eq!(b0, vec![0.0; 4]);
    let (l1, b1) = eval_bspline_basis(&[0.0, 1.0], 0, 0.5);
    assert_eq!(l1, 0);
    assert!(b1.is_empty());

    // 3. Unbounded BBC span / NaN bin_sec must not OOM/panic
    let crazy_rec = BbcRecord {
        time_ins: 1.0e15,
        lld: 1,
        jan: 1,
        ldu: 1,
        ige: 1,
        total_rate_counts: 3,
    };
    let s_crazy = BbcSeries::from_records(&[crazy_rec], 0.0, 32.768);
    assert!(s_crazy.bin_times_ins.is_empty());
    let s_nan = BbcSeries::from_records(&[crazy_rec], 0.0, f64::NAN);
    assert!(s_nan.bin_times_ins.is_empty());

    // 4. Corrupt MKF huge span in write_orbit_from_mkf must return Err (no OOM)
    let dir = tempdir().unwrap();
    let mut cols = HashMap::new();
    for k in ["POSX", "POSY", "POSZ", "VELX", "VELY", "VELZ"] {
        cols.insert(k.to_string(), vec![7000.0, 7000.0]);
    }
    let bad_mkf = MkfTable::from_columns(vec![1.0e8, 2.0e8], cols).unwrap();
    assert!(write_orbit_from_mkf(&dir.path().join("bad.orb"), &bad_mkf, &OgipMeta::default()).is_err());

    // 5. interp_lon_deg guards 359 -> 1 deg wrap discontinuity (does not return ~180 deg)
    let lon_wrap = interp_lon_deg(&[0.0, 10.0], &[359.0, 1.0], 5.0);
    assert!((lon_wrap - 359.0).abs() < 1e-9);

    // 6. Empty auto-GTI intersected with non-empty user_gti must remain empty (∅ ∩ U = ∅)
    let tct = LaxpcTct::new(
        TctTable::from_pairs(
            vec![189_900_000.0, 189_900_100.02],
            vec![1000.0, 1100.0],
        )
        .unwrap(),
    )
    .unwrap();
    let mut mkf_cols = HashMap::new();
    mkf_cols.insert("ELV".into(), vec![-20.0, -20.0]); // 100% occulted -> empty auto-GTI
    mkf_cols.insert("EARTHLON".into(), vec![100.0, 100.0]);
    let mkf_occ = MkfTable::from_columns(vec![189_900_000.0, 189_900_100.0], mkf_cols).unwrap();
    let raw_occ = RawBinnedCounts {
        t0_ins: 1000.0,
        tbin: 1.0,
        selected_counts: vec![500.0; 20],
        all_main_counts: vec![500.0; 20],
    };
    let user_gti = GtiSet::new(vec![GtiInterval::new(189_900_000.0, 189_900_050.0)]);
    let gti_empty = generate_standard_gti(
        LaxpcUnit::Lxp10,
        &raw_occ,
        &tct,
        &mkf_occ,
        None,
        &[],
        &GtiConfig::default(),
        Some(&user_gti),
    )
    .unwrap();
    assert!(gti_empty.gti_utc.is_empty());
    assert!(gti_empty.gti_ins.is_empty());

    // 7. LAXPC20 paired vs unpaired > 550 / > 2977 gain drift (`backshiftv3.f:298-307`)
    let gt = GainTable {
        entries: vec![
            GainEntry {
                day: 100.0,
                peak_30kev: [240.0, 224.0, 245.0],
                fwhm_30kev: [15.0; 3],
                peak_60kev: [470.0; 3],
                fwhm_60kev: [15.0; 3],
            },
            GainEntry {
                day: 600.0,
                peak_30kev: [235.0, 220.0, 240.0],
                fwhm_30kev: [15.0; 3],
                peak_60kev: [460.0; 3],
                fwhm_60kev: [15.0; 3],
            },
            GainEntry {
                day: 3000.0,
                peak_30kev: [230.0, 215.0, 0.0],
                fwhm_30kev: [15.0; 3],
                peak_60kev: [450.0; 3],
                fwhm_60kev: [15.0; 3],
            },
        ],
    };
    // Single-epoch query uses gaina8.dat verbatim (220.0, NOT 230.1 - 0.00571*600 = 226.674)
    assert!((gt.peak_30kev_at(LaxpcUnit::Lxp20, 600.0, None).unwrap() - 220.0).abs() < 1e-9);
    // Unpaired (source=600, bkg=100) uses gaina8.dat (220.0 vs 224.0)
    let sol_unpaired =
        compute_gain_shift(LaxpcUnit::Lxp20, 600.0, 100.0, None, &gt, None, false, false).unwrap();
    assert!((sol_unpaired.source_peak_ch - 220.0).abs() < 1e-9);
    assert!((sol_unpaired.background_peak_ch - 224.0).abs() < 1e-9);
    // Paired (> 550) uses linear drift fit (230.1 - 0.00571*day)
    let sol_paired =
        compute_gain_shift(LaxpcUnit::Lxp20, 700.0, 600.0, None, &gt, None, false, false).unwrap();
    assert!((sol_paired.source_peak_ch - (230.1 - 0.00571 * 700.0)).abs() < 1e-9);
    assert!((sol_paired.background_peak_ch - (230.1 - 0.00571 * 600.0)).abs() < 1e-9);
    // Paired (> 2977) uses post-HV-reset linear drift fit
    let sol_post =
        compute_gain_shift(LaxpcUnit::Lxp20, 3050.0, 3000.0, None, &gt, None, false, false)
            .unwrap();
    assert!((sol_post.source_peak_ch - (228.2 - 0.00571 * (3050.0 - 2977.0))).abs() < 1e-9);
    assert!((sol_post.background_peak_ch - (228.2 - 0.00571 * (3000.0 - 2977.0))).abs() < 1e-9);

    // 8. BackgroundCatalog resolves exact HV-jump background peaks (`#81` aug17 LX30 & apr18 LX10)
    let ref_dir = std::env::var("ASTROSAT_TEST_DATA")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|_| std::path::PathBuf::from("/mnt/d/ISRO/data"))
        .join("ref_src/laxpc");
    let cat_builtin = BackgroundCatalog::builtin();
    assert_eq!(
        cat_builtin.lookup_peak_30kev(LaxpcUnit::Lxp30, Some(Path::new("backlxp30aug17.pha")), 694.0),
        Some(213.0)
    );
    assert_eq!(
        cat_builtin.lookup_peak_30kev(LaxpcUnit::Lxp10, Some(Path::new("backlxp10apr18.pha")), 917.0),
        Some(157.8)
    );
    // Day-only fallback (renamed file like `lxp3level2back.pha`) still matches `|ndayb - ndy| < 2`
    assert_eq!(
        cat_builtin.lookup_peak_30kev(LaxpcUnit::Lxp30, Some(Path::new("lxp3level2back.pha")), 694.0),
        Some(213.0)
    );
    if ref_dir.is_dir() {
        let cat_file = BackgroundCatalog::load_or_builtin(&[Some(ref_dir.as_path())]);
        let gt_real = GainTable::from_file(&ref_dir.join("gaina8.dat")).unwrap();
        let hv_real = HvJumpTable::from_file(&ref_dir.join("lxhv.dat")).unwrap();

        let p_aug17 = ref_dir.join("backlxp30aug17.pha");
        let chb_aug17 = cat_file.lookup_peak_30kev(LaxpcUnit::Lxp30, Some(&p_aug17), 694.0);
        assert_eq!(chb_aug17, Some(213.0));
        let sol_aug17 = compute_gain_shift(
            LaxpcUnit::Lxp30,
            700.0,
            694.0,
            chb_aug17,
            &gt_real,
            Some(&hv_real),
            false,
            false,
        )
        .unwrap();
        assert!(
            (sol_aug17.delta_ch_30kev - 29.8).abs() < 1.0,
            "expected delta_ch ~ +29.8 for LX30 aug17->d700, got {}",
            sol_aug17.delta_ch_30kev
        );

        let p_apr18 = ref_dir.join("backlxp10apr18.pha");
        let chb_apr18 = cat_file.lookup_peak_30kev(LaxpcUnit::Lxp10, Some(&p_apr18), 917.0);
        assert_eq!(chb_apr18, Some(157.8));
        let sol_apr18 = compute_gain_shift(
            LaxpcUnit::Lxp10,
            950.0,
            917.0,
            chb_apr18,
            &gt_real,
            Some(&hv_real),
            false,
            false,
        )
        .unwrap();
        assert!(
            (sol_apr18.delta_ch_30kev - (-55.3)).abs() < 2.0,
            "expected delta_ch ~ -55.3 for LX10 apr18->d950, got {}",
            sol_apr18.delta_ch_30kev
        );

        // 9. Post-day-917 LAXPC10 shifted .channel map selection (`#105` backshiftv3.f:334-346)
        assert_eq!(ChannelMap::filename_for_day(LaxpcUnit::Lxp10, 900.0), "lxp10.channel");
        assert_eq!(ChannelMap::filename_for_day(LaxpcUnit::Lxp10, 950.0), "lxp10cshm141v2.2.channel");
        assert_eq!(ChannelMap::filename_for_day(LaxpcUnit::Lxp10, 1000.0), "lxp10cshm175v2.2.channel");
        assert_eq!(ChannelMap::filename_for_day(LaxpcUnit::Lxp10, 1050.0), "lxp10cshm197v2.2.channel");
        assert_eq!(ChannelMap::filename_for_day(LaxpcUnit::Lxp10, 1200.0), "lxp10cshm212v2.2.channel");
        let cmap_950 = ChannelMap::from_caldb_for_day(&ref_dir, LaxpcUnit::Lxp10, 950.0).unwrap();
        let pha_apr18 = BackgroundPha::from_file(&p_apr18).unwrap();
        let sh_apr18 = shift_background_spectrum(
            LaxpcUnit::Lxp10,
            sol_apr18.delta_ch_30kev,
            sol_apr18.secular_scale,
            950.0,
            &cmap_950,
            &pha_apr18.rates,
        )
        .unwrap();
        let (_cdel, fshbk) = background_lightcurve_shift_fraction(
            &pha_apr18.rates,
            &sh_apr18,
            &AnodeSelection::AllMain,
            0,
            511,
        );
        assert!(
            (fshbk - 0.550028).abs() < 1e-4,
            "expected fshbk ~ 0.550028 with lxp10cshm141v2.2.channel, got {fshbk}"
        );
    }
}
