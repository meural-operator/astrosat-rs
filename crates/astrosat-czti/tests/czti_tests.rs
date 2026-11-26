//! Unit and real-data integration tests for `astrosat-czti`.

use std::path::Path;

use astrosat_czti::bindata::*;
use astrosat_czti::bunchclean::*;
use astrosat_czti::datasel::*;
use astrosat_czti::evtclean::*;
use astrosat_czti::evtencor::*;
use astrosat_czti::fits_io::*;
use astrosat_czti::gaas::*;
use astrosat_czti::gtigen::*;
use astrosat_czti::image::*;
use astrosat_czti::noiseclean::*;
use astrosat_czti::pha2energy::*;
use astrosat_czti::pipeline::*;
use astrosat_czti::pixclean::*;
use astrosat_czti::science2event::*;
use astrosat_czti::types::*;
use astrosat_fits::{validate_ogip, FitsFile};
use tempfile::tempdir;

/// Helper: encode 1024 words into a 2048-byte big-endian CZTI packet.
fn words_to_be_packet(words: &[u16; PACKET_WORDS]) -> [u8; PACKET_SIZE_BYTES] {
    let mut out = [0u8; PACKET_SIZE_BYTES];
    for (i, &w) in words.iter().enumerate() {
        let b = w.to_be_bytes();
        out[2 * i] = b[0];
        out[2 * i + 1] = b[1];
    }
    out
}

#[test]
fn detxy_mapping_and_inverse_all_quadrants() {
    for q in 0..4u8 {
        for det in 0..16u8 {
            for pix in 0..=255u8 {
                let (detx, dety) = detid_pixid_to_detxy(q, det, pix);
                assert!(detx < 64 && dety < 64);
                let (d2, p2) = detxy_to_detid_pixid(q, detx, dety);
                assert_eq!((d2, p2), (det, pix), "roundtrip failed for q={q}, det={det}, pix={pix}");

                let (lx, ly) = detxy_to_locxy(q, detx, dety);
                assert!(lx < 128 && ly < 128);
            }
        }
    }
    // Check explicit flip rules: Q0 flips dety, Q1 flips detx
    assert_eq!(detid_pixid_to_detxy(0, 0, 0), (0, 63));
    assert_eq!(detid_pixid_to_detxy(3, 0, 0), (0, 63));
    assert_eq!(detid_pixid_to_detxy(1, 0, 0), (63, 0));
    assert_eq!(detid_pixid_to_detxy(2, 0, 0), (63, 0));
}

#[test]
fn science2event_normal_reduced_and_onboard_bunch_decoding() {
    // 1. Construct a Mode 0 (3-word normal event + VetoSpec) packet 0 with 2 normal events + 1 onboard bunch
    let mut w = [0u16; PACKET_WORDS];
    w[0] = SYNC1;
    w[1] = SYNC2;
    // mml=0, dataID=0 (Q0), pktNo=0, modeID=0
    w[2] = 0x0000;
    // 12 frame header words + 240 feb+veto + 5*3 event words (2 normal + 3-record onboard bunch)
    w[3] = 12 + 240 + 15;
    w[4] = 42; // frameNo
    w[10] = 0x0001; // secondcount hi
    w[11] = 0x86A0; // secondcount lo -> 0x000186A0 = 100000 s
    w[15] = (2 << 14) | 100; // bootpageno = 2 (onboard bunch enabled), dcnt = 100

    // FEB header at word 16: EventReadMode=1 (bit 5), VetoSpecRange=0
    w[16] = 1 << 5;
    // Veto spectrum starts at word 16 + 8 = 24
    w[24] = 777;

    // Event records start at word 16 + 240 = 256
    // Event 1: normal, ntick=500 (0.01 s), pha=321, detid=5, pixid=42, veto=3, alpha=0
    w[256] = 500;
    w[257] = ((321 << 2) << 4) | 5;
    w[258] = (42 << 8) | (3 << 1);

    // Event 2: normal, ntick=1000 (0.02 s), pha=600, detid=7, pixid=100, veto=0, alpha=1
    w[259] = 1000;
    w[260] = ((600 << 2) << 4) | 7;
    w[261] = (100 << 8) | 1;

    // Event 3..5: onboard bunch (bunch_sync=7 in top 3 bits of word 0)
    // num_bunchevents = 4 (5 events), time_dsl = 3 (60 us), time_dfs = 1 (20 us)
    w[262] = (7 << 13) | (4 << 7) | (3 << 1) | 1;
    w[263] = ((150 << 2) << 4) | 2; // detid=2, pha=150
    w[264] = 11 << 8; // pixid=11
    // Second record of bunch: ntick=2000
    w[265] = 2000;
    w[266] = ((160 << 2) << 4) | 2;
    w[267] = 12 << 8;
    // Third record of bunch: detid1..4 = (2, 3, 2, 4)
    w[268] = 2 | (3 << 4) | (2 << 8) | (4 << 12);
    w[269] = ((170 << 2) << 4) | 3;
    w[270] = 13 << 8;

    let raw = words_to_be_packet(&w);
    let stream = decode_quad_packets(0, &raw, &Science2EventConfig::default(), |t| t + 1000.0);

    assert_eq!(stream.total_packets, 1);
    assert_eq!(stream.valid_packets, 1);
    assert_eq!(stream.events.len(), 5);
    assert_eq!(stream.events.pha[0], 321);
    assert_eq!(stream.events.detid[0], 5);
    assert_eq!(stream.events.pixid[0], 42);
    assert_eq!(stream.events.veto[0], 3);
    assert_eq!(stream.events.alpha[0], 0);
    assert!((stream.events.cztseccnt[0] - 100_000.01).abs() < 1e-9);
    assert!((stream.events.time[0] - 101_000.01).abs() < 1e-9);

    assert_eq!(stream.events.pha[1], 600);
    assert_eq!(stream.events.alpha[1], 1);

    // Onboard bunch checks
    assert_eq!(stream.bunches.len(), 1);
    assert_eq!(stream.bunches.row_num[0], 2);
    assert_eq!(stream.bunches.num_event[0], 4);
    assert_eq!(stream.bunches.time_dfs[0], 1);
    assert_eq!(stream.bunches.time_dsl[0], 3);
    assert_eq!(stream.bunches.detid1[0], 2);
    assert_eq!(stream.bunches.detid2[0], 3);
    assert_eq!(stream.bunches.detid4[0], 4);

    // VetoSpec check (vetoSpecRange == 0 -> offset 24)
    assert_eq!(stream.veto_spec.veto_spec.len(), 1);
    assert_eq!(stream.veto_spec.veto_spec[0][24], 777);

    // 2. Test reduced mode (modeID=3, 2-word events, feb_veto_size=8)
    let mut wr = [0u16; PACKET_WORDS];
    wr[0] = SYNC1;
    wr[1] = SYNC2;
    wr[2] = (1 << 8) | 3; // dataID=1 (Q1), pktNo=0, modeID=3
    wr[3] = 12 + 8 + 2; // 1 reduced event
    wr[11] = 500;
    wr[16] = 1 << 5; // EventReadMode=1
    // Reduced event at word 16 + 8 = 24: ntick=100 (9 bits << 7), pha_hi=25 (7 bits)
    // word 25: pha_lo=2 (2 bits << 14), pixid=88 (8 bits << 6), veto=1 (bit 5), alpha=0 (bit 4), detid=9
    wr[24] = (100 << 7) | 25;
    wr[25] = (2 << 14) | (88 << 6) | (1 << 5) | 9;
    let raw_r = words_to_be_packet(&wr);
    let stream_r = decode_quad_packets(1, &raw_r, &Science2EventConfig::default(), |t| t);
    assert_eq!(stream_r.events.len(), 1);
    assert_eq!(stream_r.events.cztntick[0], 100);
    assert_eq!(stream_r.events.pha[0], ((25 << 2) | 2) << 1);
    assert_eq!(stream_r.events.detid[0], 9);
    assert_eq!(stream_r.events.pixid[0], 88);
    assert_eq!(stream_r.events.veto[0], 1);
}

#[test]
fn bunchclean_removes_bunches_100s_ticks_and_post_bunch_deadtime() {
    let mut ev = QuadEvents::new(0);
    let push = |ev: &mut QuadEvents, t: f64, czt: f64, det: u8, pix: u8| {
        let (dx, dy) = detid_pixid_to_detxy(0, det, pix);
        ev.time.push(t);
        ev.cztseccnt.push(czt);
        ev.cztntick.push(0);
        ev.pha.push(100);
        ev.detid.push(det);
        ev.pixid.push(pix);
        ev.detx.push(dx);
        ev.dety.push(dy);
        ev.veto.push(0);
        ev.alpha.push(0);
    };

    // Event 0: inside first 400 ms of 100s CZT clock tick (100.10 s) -> MUST be rejected!
    push(&mut ev, 1000.10, 100.10, 1, 10);
    // Event 1: clean single at 1001.00 s -> kept
    push(&mut ev, 1001.00, 101.00, 1, 10);
    // Events 2,3: clean double (20 us apart) at 1002.00 s -> both kept
    push(&mut ev, 1_002.000_000, 102.000_000, 2, 20);
    push(&mut ev, 1_002.000_020, 102.000_020, 2, 21);
    // Events 4..7: 4-event cosmic-ray bunch (20 us apart, dominant det=3) -> all removed
    push(&mut ev, 1_003.000_000, 103.000_000, 3, 1);
    push(&mut ev, 1_003.000_020, 103.000_020, 3, 2);
    push(&mut ev, 1_003.000_040, 103.000_040, 3, 3);
    push(&mut ev, 1_003.000_060, 103.000_060, 4, 4);
    // Event 8: 0.5 ms after bunch in SAME dominant det=3 (< skipT2=1ms) -> rejected!
    push(&mut ev, 1_003.000_560, 103.000_560, 3, 5);
    // Event 9: 0.5 ms after bunch in DIFFERENT det=5 (skipT1=0) -> kept!
    push(&mut ev, 1_003.000_600, 103.000_600, 5, 6);
    // Event 10: sentinel event at 1005.0 s (loop stops at nrows-1 in cztbunchcleanv2)
    push(&mut ev, 1_005.000_000, 105.000_000, 6, 7);

    let res = bunchclean_quad(
        &ev,
        &BunchTable::new(0),
        &[(1000.0, 1006.0)],
        &BunchCleanConfig::default(),
    );
    assert_eq!(res.out_bunches.len(), 1);
    assert_eq!(res.out_bunches.num_event[0], 3); // 4 - 1
    assert_eq!(res.events.len(), 4); // ev 1, 2, 3, 9
    assert_eq!(res.events.detid, vec![1, 2, 2, 5]);
    // Dominant module 3 suffered 1 ms deadtime + 60 us bunch -> lower exposure than module 5
    assert!(res.pixel_exposure[3 * 256] < res.pixel_exposure[5 * 256]);
}

#[test]
fn pha2energy_5_temp_bins_and_ebounds() {
    let mut ev = QuadEvents::new(0);
    // 5 events at different times spanning temperatures 0, 5, 10, 15, 20 C, plus PHA==0
    for (i, t) in [0.0, 10.0, 20.0, 30.0, 40.0, 50.0].iter().enumerate() {
        ev.time.push(*t);
        ev.cztseccnt.push(*t);
        ev.cztntick.push(0);
        ev.pha.push(if i == 5 { 0 } else { 100 });
        ev.detid.push(0);
        ev.pixid.push(0);
        ev.detx.push(0);
        ev.dety.push(63);
        ev.veto.push(0);
        ev.alpha.push(0);
    }

    let mut temp_tbl = TempTable::default();
    for (t, tc) in [(0.0, 1.0f32), (10.0, 5.0), (20.0, 10.0), (30.0, 15.0), (40.0, 20.0)] {
        temp_tbl.time.push(t);
        temp_tbl.cztseccnt.push(t);
        temp_tbl.temperature.push([tc; 16]);
        temp_tbl.quad_flag.push(0);
    }
    let tseries = QuadTempSeries::from_temp_table(&temp_tbl, 0);

    let mut go = QuadGainOffset::default();
    for b in 0..5 {
        go.gain[b][0] = 0.2 + 0.1 * (b as f32);
        go.offset[b][0] = 5.0;
    }
    let ebounds = EboundsTable::default();

    pha2energy_quad(&mut ev, &tseries, &go, &ebounds);
    assert!((ev.energy[0] - 25.0).abs() < 1e-4);
    assert!((ev.energy[1] - 35.0).abs() < 1e-4);
    assert!((ev.energy[2] - 45.0).abs() < 1e-4);
    assert!((ev.energy[3] - 55.0).abs() < 1e-4);
    assert!((ev.energy[4] - 65.0).abs() < 1e-4);
    assert_eq!(ev.energy[5], 0.0);
    assert_eq!(ev.pi[0], 40); // (25 - 5) / 0.5 = 40
    assert_eq!(ev.pi[5], 0);
}

#[test]
fn gtigen_gaas_datasel_pixclean_noiseclean_evtencor_evtclean() {
    // 1. gtigen
    let rules = parse_mkf_threshold_text("ELV 1 5 -\nANG_OFFSET 1 - 0.03\nMOON_ANGLE 0 - -\n");
    assert_eq!(rules.len(), 2);
    let gtis = generate_czti_gtis(
        &[(0.0, 100.0)],
        &[(10.0, 10.5), (20.0, 80.0)],
        &[
            vec![(20.0, 80.0)],
            vec![(20.0, 80.0)],
            vec![(20.0, 80.0)],
            vec![(20.0, 80.0)],
        ],
        Some(&[(40.0, 45.0)]),
        Some(&[(15.0, 70.0)]),
        1.05,
    )
    .unwrap();
    // [10, 10.5] dropped (< 1.05s); [20, 80] - [40, 45] & [15, 70] -> [20, 40], [45, 70]
    assert_eq!(gtis.common, vec![(20.0, 40.0), (45.0, 70.0)]);

    // 2. gaas: identity alignment + satellite yaw=X, roll=Y, pitch=Z -> pointing (0,0,1) => pitch axis (0,0,1) => DEC=90
    let rpy = vec![RpySampleRad {
        time: 20.0,
        yaw_ra: 0.0,
        yaw_dec: 0.0,
        roll_ra: std::f64::consts::FRAC_PI_2,
        roll_dec: 0.0,
        pitch_ra: 0.0,
        pitch_dec: std::f64::consts::FRAC_PI_2,
    }];
    let asp = compute_quad_aspect(0, "test.teldef".into(), &TeldefMatrix::default(), &rpy, 20.0, 25.0, None);
    assert!((asp.dec_deg - 90.0).abs() < 1e-3, "got dec={}", asp.dec_deg);

    // 3. evtencor secular gain drift & idempotency
    let drift_at_epoch = secular_gain_drift_factor(GAIN_DRIFT_EPOCH_SEC, true);
    assert!((drift_at_epoch - 0.994).abs() < 1e-6);
    let drift_1yr = secular_gain_drift_factor(GAIN_DRIFT_EPOCH_SEC + 365.0 * 86400.0, true);
    assert!(drift_1yr > drift_at_epoch);

    let mut ev = QuadEvents::new(0);
    for i in 0..10 {
        ev.time.push(25.0 + (i as f64));
        ev.cztseccnt.push(25.0 + (i as f64));
        ev.cztntick.push(0);
        ev.pha.push(200);
        ev.detid.push(0);
        ev.pixid.push(i as u8);
        ev.detx.push(i as u8);
        ev.dety.push(63);
        ev.veto.push(if i == 8 { 5 } else { 0 });
        ev.alpha.push(if i == 9 { 1 } else { 0 });
        ev.pi.push(100);
        ev.energy.push(55.0);
    }
    let ds = datasel_quad(&ev, &gtis.common);
    assert_eq!(ds.len(), 10);

    let mut gc = QuadGainCorCaldb::default();
    gc.lld[0] = 120; // pixel 0 should be cut by LLD!
    let ec1 = evtencor_quad(&ds, &gc, GAIN_DRIFT_EPOCH_SEC, false, &EvtEnCorConfig::default());
    assert_eq!(ec1.removed_by_lld_uld, 1);
    assert_eq!(ec1.events.len(), 9);
    // Idempotency: running again with already_gain_corrected=true does NOTHING
    let ec2 = evtencor_quad(&ec1.events, &gc, GAIN_DRIFT_EPOCH_SEC, true, &EvtEnCorConfig::default());
    assert!(ec2.skipped_already_corrected);
    assert_eq!(ec2.events.pi, ec1.events.pi);

    // 4. evtclean drops alpha != 0 (i=9) and veto != 0 (i=8)
    let cleaned = evtclean_quad(&ec1.events, &EvtCleanConfig::default());
    assert_eq!(cleaned.len(), 7);

    // 5. eventsep separates a 35 us double from singles
    let mut ev_sep = QuadEvents::new(0);
    for (i, &t) in [1.0, 2.0, 2.000_035, 3.0].iter().enumerate() {
        ev_sep.time.push(t);
        ev_sep.cztseccnt.push(t);
        ev_sep.cztntick.push(0);
        ev_sep.pha.push(100);
        ev_sep.detid.push(0);
        ev_sep.pixid.push(i as u8);
        ev_sep.detx.push(i as u8);
        ev_sep.dety.push(63);
        ev_sep.veto.push(0);
        ev_sep.alpha.push(0);
    }
    let (singles, doubles) = eventsep_quad(&ev_sep, 40e-6);
    assert_eq!(singles.len(), 2);
    assert_eq!(doubles.len(), 2);

    // 6. pixclean 5-sigma iterative clipping flags a hot pixel
    let mut dph = [10u64; PIXELS_PER_QUAD];
    dph[42] = 5000;
    dph[99] = 0;
    let mut flags = [pix_flag::GOOD; PIXELS_PER_QUAD];
    find_bad_pixels(&dph, &mut flags, 5.0);
    assert_eq!(flags[42], pix_flag::NOISY);
    assert_eq!(flags[99], pix_flag::DEAD);
    assert_eq!(flags[0], pix_flag::GOOD);
}

#[test]
fn radec2czti_mask_weight_normalization_imaging_and_ogip_outputs() {
    // 1. On-axis source at Roll=(83.63, 22.01, rot=0) -> tx=0, ty=0
    let (tx0, ty0) = radec_to_czti_tx_ty(83.63, 22.01, 0.0, 83.63, 22.01);
    assert!(tx0.abs() < 1e-6 && ty0.abs() < 1e-6, "got ({tx0}, {ty0})");
    assert!(validate_fov(tx0, ty0).is_ok());

    // Off-axis source 10 deg away -> rejected by FOV check (> 5.8 deg)
    let (tx_far, ty_far) = radec_to_czti_tx_ty(83.63, 22.01, 0.0, 94.0, 22.01);
    assert!(validate_fov(tx_far, ty_far).is_err());

    // 2. Coded-mask weight renormalization: uniform background must cancel to ZERO!
    let mut mask64 = vec![vec![0.0f32; QUAD_DIM]; QUAD_DIM];
    let mut rng_state: u64 = 0x5DEECE66D1234567;
    for row in mask64.iter_mut() {
        for cell in row.iter_mut() {
            // Deterministic aperiodic xorshift64 coded mask (~50% open)
            rng_state ^= rng_state << 13;
            rng_state ^= rng_state >> 7;
            rng_state ^= rng_state << 17;
            *cell = if (rng_state & 1) == 0 { 1.0 } else { 0.0 };
        }
    }
    let (open_frac, raw_w) = raytrace_quad_shadow_simple(0.0, 0.0, &mask64);
    let bp = BadpixTable::default();
    let bkgscal = vec![1.0f32; PIXELS_PER_QUAD * PI_CHANNELS];
    let lld_map = [[0i16; QUAD_DIM]; QUAD_DIM];
    let uld_map = [[511i16; QUAD_DIM]; QUAD_DIM];
    let pix_exp = [1.0f64; PIXELS_PER_QUAD];
    let mw = calculate_renormalized_weights(
        0, &open_frac, &raw_w, &bp, 0, &bkgscal, &lld_map, &uld_map, &pix_exp, 100,
    );
    // Sum of renormalized weights * BKGSCAL * pix_exp across a module MUST be ~0 (exact background cancellation!)
    let mut bkg_sum = 0.0f64;
    for dety in 0..16 {
        for detx in 0..16 {
            bkg_sum += mw.weight(detx, dety, 100) as f64;
        }
    }
    assert!(
        bkg_sum.abs() < 1e-4,
        "renormalized weights failed background cancellation: sum={bkg_sum}"
    );

    // 3. Sky imaging (`balcc` and `fft`) recovers on-axis source at center
    let mut dpi = vec![vec![0.0f32; QUAD_DIM]; QUAD_DIM];
    for y in 0..QUAD_DIM {
        for x in 0..QUAD_DIM {
            dpi[y][x] = 10.0 + 50.0 * mask64[y][x];
        }
    }
    let asp = QuadAspect {
        quad_id: 0,
        ra_deg: 83.63,
        dec_deg: 22.01,
        twist_deg: 0.0,
        ..Default::default()
    };
    let img_balcc = balanced_cross_correlation_from_mask(&dpi, &mask64, 5, &asp);
    let (px, py, pv) = img_balcc.peak().unwrap();
    assert_eq!((px, py), (5, 5), "balcc peak should be centered at (5,5), got ({px},{py})");
    assert!((pv - 50.0).abs() < 1.0, "balcc amplitude expected ~50, got {pv}");

    let img_fft = fft_cross_correlation(&dpi, &mask64, 1, &asp);
    let (fx, fy, _) = img_fft.peak().unwrap();
    assert_eq!((fx, fy), (img_fft.width / 2, img_fft.height / 2));

    // 3b. Off-axis FFT sign convention & Q2/Q3 yshift axis verification:
    // A source at (+tx, -ty) casts a mask shadow shifted on the detector by (-dx, +dy)
    // (i.e. `dpi[y][x] = mask[y - dy][x + dx]` with dx=4, dy=3).
    // Legacy `fft.h` places the sky peak at `(center + dx, center - dy)` (NOT 180-deg inverted!).
    let mut dpi_off = vec![vec![0.0f32; QUAD_DIM]; QUAD_DIM];
    for (y, row) in dpi_off.iter_mut().enumerate() {
        for (x, cell) in row.iter_mut().enumerate() {
            let sx = x as isize + 4;
            let sy = y as isize - 3;
            let m = if sx >= 0 && sx < QUAD_DIM as isize && sy >= 0 && sy < QUAD_DIM as isize {
                mask64[sy as usize][sx as usize]
            } else {
                0.0
            };
            *cell = 10.0 + 50.0 * m;
        }
    }
    let img_fft_off = fft_cross_correlation_with_shift(&dpi_off, &mask64, 1, 0.0, 0.0, &asp);
    let (ofx, ofy, _) = img_fft_off.peak().unwrap();
    let cx = img_fft_off.width / 2;
    let cy = img_fft_off.height / 2;
    assert_eq!(
        (ofx, ofy),
        (cx + 4, cy - 3),
        "FFT off-axis source sign should be (+4, -3), got ({}, {})",
        ofx as isize - cx as isize,
        ofy as isize - cy as isize
    );

    // Passing `yshift = 1.0` (as for Q2/Q3) must shift the peak along Y (not X!) by -1*ov.
    let img_fft_ysh = fft_cross_correlation_with_shift(&dpi, &mask64, 2, 0.0, 1.0, &asp);
    let (yx, yy, _) = img_fft_ysh.peak().unwrap();
    assert_eq!(
        (yx, yy),
        (img_fft_ysh.width / 2, img_fft_ysh.height / 2 - 2),
        "FFT yshift=1.0 (ov=2) should shift Y by -2 and leave X centered, got ({yx},{yy})"
    );

    // 3c. SHADOWLIB `balanced_cross_correlation` layout (ThetaX outer, ThetaY inner) + 0.10 deg WCS
    let grid_dim = 61usize;
    let target_tx_i = 54usize; // +2.40 deg
    let target_ty_i = 12usize; // -1.80 deg
    let mut synth_shadows: Vec<Vec<f32>> = vec![vec![0.5f32; PIXELS_PER_QUAD]; grid_dim * grid_dim];
    let mut flat_mask = vec![0.0f32; PIXELS_PER_QUAD];
    for (y, row) in mask64.iter().enumerate() {
        for (x, &v) in row.iter().enumerate() {
            flat_mask[y * QUAD_DIM + x] = v;
        }
    }
    // SHADOWLIB row index = tx_i * grid_dim + ty_i
    synth_shadows[target_tx_i * grid_dim + target_ty_i] = flat_mask;
    let shadow_refs: Vec<&[f32]> = synth_shadows.iter().map(|v| v.as_slice()).collect();
    let img_sl = balanced_cross_correlation(&dpi, &shadow_refs, grid_dim, &asp);
    let (slx, sly, _) = img_sl.peak().unwrap();
    assert_eq!(
        (slx, sly),
        (target_tx_i, target_ty_i),
        "SHADOWLIB balcc peak should be at (tx_i={target_tx_i}, ty_i={target_ty_i}), got ({slx},{sly})"
    );
    assert!(
        (img_sl.wcs.cdelt1 - 0.10).abs() < 1e-12 && (img_sl.wcs.cdelt2 - 0.10).abs() < 1e-12,
        "SHADOWLIB balcc WCS scale should be 0.10 deg/pix, got ({}, {})",
        img_sl.wcs.cdelt1,
        img_sl.wcs.cdelt2
    );
    let tx_wcs = (slx as f64 + 1.0 - img_sl.wcs.crpix1b) * img_sl.wcs.cdelt1b;
    let ty_wcs = (sly as f64 + 1.0 - img_sl.wcs.crpix2b) * img_sl.wcs.cdelt2b;
    assert!(
        (tx_wcs - 2.40).abs() < 1e-9 && (ty_wcs + 1.80).abs() < 1e-9,
        "SHADOWLIB WCS (TX, TY) expected (+2.40, -1.80), got ({tx_wcs}, {ty_wcs})"
    );

    // 3d. Legacy `.pha` QUALITY & SNR-based logarithmic GROUPING parity (`level2validation.cpp`)
    let mut test_cnt = vec![0.0f32; PI_CHANNELS];
    let mut test_err = vec![0.0f32; PI_CHANNELS];
    for ch in 10..=450 {
        test_cnt[ch] = 100.0;
        test_err[ch] = 10.0;
    }
    let (q_flags, grp_flags) = compute_spectrum_quality_and_grouping(&test_cnt, &test_err, true);
    assert!(q_flags[0..34].iter().all(|&q| q == 5), "channels < 34 must have QUALITY=5");
    assert!(q_flags[34..=390].iter().all(|&q| q == 0), "good channels 34..=390 must have QUALITY=0");
    assert!(q_flags[391..PI_CHANNELS].iter().all(|&q| q == 5), "channels > 390 must have QUALITY=5");
    assert!(grp_flags[0..34].iter().all(|&g| g == 1), "ungrouped low channels must be 1");
    assert!(grp_flags[34..=390].contains(&-1), "SNR-based logarithmic grouping should group high-E bins");
    assert!(grp_flags[391..PI_CHANNELS].iter().all(|&g| g == 1), "ungrouped high channels must be 1");
    let (q_zero, grp_zero) =
        compute_spectrum_quality_and_grouping(&vec![0.0; PI_CHANNELS], &vec![0.0; PI_CHANNELS], true);
    assert!(q_zero.iter().all(|&q| q == 0), "zero-count spectrum keeps QUALITY=0");
    assert!(grp_zero.iter().all(|&g| g == 1), "zero-count spectrum keeps GROUPING=1");

    // 4. Tantalum absorption sanity check (K-edge jump at 67.42 keV)
    assert!(absco_ta(70.0) > absco_ta(60.0));
    assert!(ta_mask_transmission(30.0) < 0.01);
    assert!(ta_mask_transmission(200.0) > 0.3);

    // 5. OGIP .pha, .lc, .rsp write & validate_ogip
    let dir = tempdir().unwrap();
    let meta = CztiObsMeta {
        obs_id: Some("9000000406".into()),
        object: Some("TEST_SRC".into()),
        tstart: 197_000_000.0,
        tstop: 197_000_100.0,
        ..Default::default()
    };
    let spec = CztiSpectrum {
        quad_id: Some(0),
        channel: (0..512).collect(),
        counts: vec![10.0; 512],
        stat_err: vec![3.16; 512],
        sys_err: vec![0.02; 512],
        quality: vec![0; 512],
        grouping: vec![1; 512],
        exposure: 100.0,
        thetax_deg: 0.0,
        thetay_deg: 0.0,
        mask_weighted: true,
        badpix_threshold: 0,
        gti_type: GtiType::Quad,
        respfile: "test.rsp".into(),
        tstart: meta.tstart,
        tstop: meta.tstop,
    };
    let pha_path = dir.path().join("test.pha");
    write_spectrum_file(&pha_path, &meta, &spec).unwrap();
    let pha_fits = FitsFile::open(&pha_path).unwrap();
    let pha_issues = validate_ogip(&pha_fits);
    assert!(pha_issues.is_empty(), "PHA OGIP issues: {pha_issues:?}");

    let lc = CztiLightCurve {
        quad_id: Some(0),
        time: vec![197_000_005.0, 197_000_015.0],
        rate: vec![12.0, 11.5],
        error: vec![1.1, 1.0],
        fracexp: vec![0.95, 0.90],
        timedel: 10.0,
        tstart: meta.tstart,
        tstop: meta.tstop,
        emin_kev: 20.0,
        emax_kev: 100.0,
        thetax_deg: 0.0,
        thetay_deg: 0.0,
        mask_weighted: true,
        badpix_threshold: 0,
    };
    let lc_path = dir.path().join("test.lc");
    write_lightcurve_file(&lc_path, &meta, &lc).unwrap();
    let lc_fits = FitsFile::open(&lc_path).unwrap();
    let lc_issues = validate_ogip(&lc_fits);
    assert!(lc_issues.is_empty(), "LC OGIP issues: {lc_issues:?}");

    let ebounds = EboundsTable::default();
    let rsp = CztiResponse {
        quad_id: Some(0),
        n_ebins: 4,
        n_pibins: 512,
        energ_lo: vec![10.0, 20.0, 30.0, 40.0],
        energ_hi: vec![20.0, 30.0, 40.0, 50.0],
        matrix: {
            let mut m = vec![0.0f32; 4 * 512];
            m[20] = 100.0;
            m[512 + 40] = 120.0;
            m[2 * 512 + 60] = 110.0;
            m[3 * 512 + 80] = 90.0;
            m
        },
        ebounds_channel: ebounds.channel,
        ebounds_emin: ebounds.e_min,
        ebounds_emax: ebounds.e_max,
    };
    let rsp_path = dir.path().join("test.rsp");
    write_response_file(&rsp_path, &meta, &rsp).unwrap();
    let rsp_fits = FitsFile::open(&rsp_path).unwrap();
    let rsp_issues = validate_ogip(&rsp_fits);
    assert!(rsp_issues.is_empty(), "RSP OGIP issues: {rsp_issues:?}");
}

#[test]
fn real_obsid_9000000406_end_to_end_pipeline_and_reference_parity() {
    let data_root = std::env::var("ASTROSAT_TEST_DATA")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|_| std::path::PathBuf::from("/mnt/d/ISRO/data"));
    let base = data_root.join("sample/czti/czti_sample_data/20160331_T01_112T01_9000000406_level1");
    let ref_l2 = data_root.join("sample/czti/czti_sample_data/20160331_T01_112T01_9000000406_level2/modeM0");
    let caldb_dir = data_root.join("caldb/czti");
    if !base.exists() || !caldb_dir.exists() {
        eprintln!(
            "Skipping real ObsID 9000000406 test (sample data not present under {})",
            data_root.display()
        );
        return;
    }

    let out_dir = tempdir().unwrap();
    let mut cfg = CztiPipelineConfig::new(
        base.join("modeM0/AS1T01_112T01_9000000406cztM0_level1.fits"),
        Some(base.join("modeSS/AS1T01_112T01_9000000406cztSS_level1.fits")),
        base.join("aux/AS1T01_112T01_9000000406czt_level1.tct"),
        base.join("AS1T01_112T01_9000000406czt_level2.mkf"),
        &caldb_dir,
        out_dir.path(),
    );
    cfg.mkf_threshold_file = Some(base.join("mkfThresholds.txt"));
    // Run full 3 stages including mask-weighted spectrum, lightcurve, and response!
    cfg.bindata.time_binsize = 10.0;
    cfg.bindata.quads = vec![0, 1, 2, 3];
    cfg.image.algorithm = ImageAlgorithm::BalCC;

    let summary = run_pipeline(&cfg).expect("full CZTI L1->L2 pipeline should succeed");

    // Verify all 4 quadrants produced events, GTIs, aspects, DPH/DPI, images, spectra, LCs, RSPs
    for q in 0..4 {
        assert!(
            summary.stage1_raw_events[q] > 10_000,
            "Q{q} expected >10k raw events, got {}",
            summary.stage1_raw_events[q]
        );
        assert!(
            summary.stage1_bunchclean_events[q] > 5_000
                && summary.stage1_bunchclean_events[q] <= summary.stage1_raw_events[q],
            "Q{q} bunchclean events out of range: {}",
            summary.stage1_bunchclean_events[q]
        );
        assert!(
            summary.stage2_clean_events[q] > 1_000,
            "Q{q} clean events too low: {}",
            summary.stage2_clean_events[q]
        );
        assert!(
            summary.quad_exposures_sec[q] > 100.0,
            "Q{q} exposure too low: {}",
            summary.quad_exposures_sec[q]
        );
    }

    assert_eq!(summary.aspect_paths.len(), 4);
    assert_eq!(summary.image_paths.len(), 4);
    assert_eq!(summary.spectrum_paths.len(), 4);
    assert_eq!(summary.lightcurve_paths.len(), 4);
    assert_eq!(summary.response_paths.len(), 4);

    // Cross-check CHECKSUM/DATASUM & OGIP conformance on generated products
    for p in summary
        .spectrum_paths
        .iter()
        .chain(summary.lightcurve_paths.iter())
        .chain(summary.response_paths.iter())
    {
        let fits = FitsFile::open(p).unwrap();
        for h in fits.iter_hdus() {
            assert_eq!(h.verify_checksum(), (true, true), "bad checksum in {}", p.display());
        }
        let issues = validate_ogip(&fits);
        assert!(issues.is_empty(), "OGIP issues in {}: {issues:?}", p.display());
    }

    // Cross-check against official reference Level-2 outputs!
    if ref_l2.exists() {
        let ref_bc = FitsFile::open(ref_l2.join("AS1T01_112T01_9000000406cztM0_level2_bc.evt")).unwrap();
        let ref_clean = FitsFile::open(ref_l2.join("AS1T01_112T01_9000000406czt_level2_quad_clean.evt")).unwrap();
        for q in 0..4 {
            let ref_bc_rows = ref_bc.hdu_by_name(&format!("Q{q}")).unwrap().as_bintable().unwrap().nrows();
            let our_bc_rows = summary.stage1_bunchclean_events[q];
            let rel_bc = ((our_bc_rows as f64) - (ref_bc_rows as f64)).abs() / (ref_bc_rows.max(1) as f64);
            eprintln!("Q{q} bc.evt: ours={our_bc_rows}, ref={ref_bc_rows} (rel diff {:.3}%)", rel_bc * 100.0);
            assert!(rel_bc < 0.05, "Q{q} bc.evt row count diverged >5%: ours={our_bc_rows}, ref={ref_bc_rows}");

            let ref_clean_hdu = ref_clean.hdu_by_name(&format!("Q{q}")).unwrap();
            let ref_clean_rows = ref_clean_hdu.as_bintable().unwrap().nrows();
            let ref_exp = ref_clean_hdu.header().get_f64("EXPOSURE").unwrap_or(0.0);
            let our_clean_rows = summary.stage2_clean_events[q];
            let our_exp = summary.quad_exposures_sec[q];
            eprintln!(
                "Q{q} clean.evt: ours={our_clean_rows} (exp={our_exp:.2}s), ref={ref_clean_rows} (exp={ref_exp:.2}s)"
            );
            let rel_exp = (our_exp - ref_exp).abs() / ref_exp.max(1.0);
            assert!(rel_exp < 0.05, "Q{q} exposure diverged >5%: ours={our_exp}, ref={ref_exp}");

            // Check aspect RA/DEC/TWIST parity against official aspect_Q{q}
            // Note: official v1.0 reference averages over all 1297 s (including pre-GTI slew),
            // while our GTI-filtered GAAS averages over the 365 s on-Crab GTI (RA=83.633, DEC=22.014).
            let ref_asp_path = ref_l2.join(format!("AS1T01_112T01_9000000406cztM0_level2.aspect_Q{q}"));
            if ref_asp_path.exists() {
                let ref_asp = FitsFile::open(&ref_asp_path).unwrap();
                let h = ref_asp.hdu_by_name("ASPECT").unwrap().header();
                let ref_ra = h.get_f64("RA").unwrap();
                let ref_dec = h.get_f64("DEC").unwrap();
                let ref_tw = h.get_f64("TWIST").unwrap();
                let (our_ra, our_dec, our_tw) = summary.quad_aspect_deg[q];
                eprintln!(
                    "Q{q} aspect (GTI): ours=({our_ra:.4}, {our_dec:.4}, {our_tw:.4}), ref(all-time)=({ref_ra:.4}, {ref_dec:.4}, {ref_tw:.4})"
                );
                assert!((our_ra - ref_ra).abs() < 0.25, "Q{q} RA mismatch: ours={our_ra}, ref={ref_ra}");
                assert!((our_dec - ref_dec).abs() < 0.25, "Q{q} DEC mismatch: ours={our_dec}, ref={ref_dec}");
                assert!((our_tw - ref_tw).abs() < 0.5, "Q{q} TWIST mismatch: ours={our_tw}, ref={ref_tw}");
            }
        }
    }

    // Also verify Generalized v3.0 noise-cleaning algorithm runs cleanly through Stage 2 on real data!
    let gen_dir = tempdir().unwrap();
    let mut gen_cfg = CztiPipelineConfig::new(
        base.join("modeM0/AS1T01_112T01_9000000406cztM0_level1.fits"),
        Some(base.join("modeSS/AS1T01_112T01_9000000406cztSS_level1.fits")),
        base.join("aux/AS1T01_112T01_9000000406czt_level1.tct"),
        base.join("AS1T01_112T01_9000000406czt_level2.mkf"),
        &caldb_dir,
        gen_dir.path(),
    );
    gen_cfg.mkf_threshold_file = Some(base.join("mkfThresholds.txt"));
    gen_cfg.evt_sel_algo = EventSelectionAlgo::Generalized;
    gen_cfg.exit_stage = 2;
    let gen_sum = run_pipeline(&gen_cfg).expect("Generalized v3.0 pipeline should succeed");
    for (q, &cnt) in gen_sum.stage2_clean_events.iter().enumerate() {
        assert!(
            cnt > 1_000,
            "Generalized Q{q} clean events too low: {cnt}"
        );
    }
    eprintln!("Generalized v3.0 clean events: {:?}", gen_sum.stage2_clean_events);

    // Run NASA HEASARC fitsverify 4.22 on ALL generated Level-2 FITS products if available
    let fitsverify = Path::new("/mnt/d/ISRO/.toolchain/bin/fitsverify");
    if fitsverify.exists() {
        let mut all_products = Vec::new();
        for p in [
            &summary.bc_evt_path,
            &summary.clean_evt_path,
            &summary.gti_path,
            &summary.badpix_path,
            &summary.livetime_path,
            &summary.dph_path,
            &summary.dpi_path,
        ]
        .into_iter()
        .flatten()
        {
            all_products.push(p.clone());
        }
        all_products.extend(summary.aspect_paths.iter().cloned());
        all_products.extend(summary.image_paths.iter().cloned());
        all_products.extend(summary.spectrum_paths.iter().cloned());
        all_products.extend(summary.lightcurve_paths.iter().cloned());
        all_products.extend(summary.response_paths.iter().cloned());

        for p in &all_products {
            let out = std::process::Command::new(fitsverify)
                .arg("-e")
                .arg(p)
                .output()
                .expect("failed to invoke fitsverify");
            let stdout = String::from_utf8_lossy(&out.stdout);
            let stderr = String::from_utf8_lossy(&out.stderr);
            assert!(
                out.status.success(),
                "fitsverify failed on {}:\nSTDOUT:\n{stdout}\nSTDERR:\n{stderr}",
                p.display()
            );
        }
        eprintln!("fitsverify 4.22: 0 errors across all {} generated CZTI L2 FITS products!", all_products.len());
    }
}

#[test]
fn audit_integrity_and_legacy_parity_regression_suite() {
    // 1. Mid-bunch frame truncation must NEVER desynchronize event columns across frames (#47-1)
    let mut f1 = [0u16; PACKET_WORDS];
    f1[0] = SYNC1;
    f1[1] = SYNC2;
    f1[2] = 0x0000;
    f1[3] = 12 + 240 + 6; // 1 normal event + ONLY 1 record of a 3-record onboard bunch!
    f1[4] = 1;
    f1[11] = 100;
    f1[15] = (2 << 14) | 10;
    f1[16] = 1 << 5;
    f1[256] = 50;
    f1[257] = ((111 << 2) << 4) | 1;
    f1[258] = 10 << 8;
    // Partial bunch record 0 only:
    f1[259] = (7 << 13) | (4 << 7) | (2 << 1) | 1;
    f1[260] = ((222 << 2) << 4) | 2;
    f1[261] = 20 << 8;

    let mut f2 = [0u16; PACKET_WORDS];
    f2[0] = SYNC1;
    f2[1] = SYNC2;
    f2[2] = 0x0000;
    f2[3] = 12 + 240 + 3; // 1 normal event in frame 2
    f2[4] = 2;
    f2[11] = 101;
    f2[15] = (2 << 14) | 10;
    f2[16] = 1 << 5;
    f2[256] = 60;
    f2[257] = ((333 << 2) << 4) | 3;
    f2[258] = 30 << 8;

    let mut raw = Vec::new();
    raw.extend_from_slice(&words_to_be_packet(&f1));
    raw.extend_from_slice(&words_to_be_packet(&f2));
    let stream = decode_quad_packets(0, &raw, &Science2EventConfig::default(), |t| t);
    assert_eq!(stream.events.len(), 2);
    assert_eq!(stream.events.pha, vec![111, 333]);
    assert_eq!(stream.events.detid, vec![1, 3]);
    assert_eq!(stream.bunches.len(), 0);

    // 2. Wichmann-Hill PRNG exact legacy modulus 30319 (#55-2)
    let mut rng = WichmannHillRng::new(5, 11, 17);
    let u0 = rng.next_f64();
    let expected = (855.0 / 30269.0 + 1892.0 / 30307.0 + 2890.0 / 30319.0_f64).fract();
    assert!((u0 - expected).abs() < 1e-12, "u0={u0}, expected={expected}");

    // 3. ADDITIONAL_LLD clamping (#55-3)
    let mut lld_q1 = vec![20i16; PIXELS_PER_QUAD];
    apply_additional_lld(1, &mut lld_q1);
    assert_eq!(lld_q1[0], 20);
    assert_eq!(lld_q1[7 * 256], 58);
    assert_eq!(lld_q1[11 * 256 + 10], 62);
    assert_eq!(lld_q1[15 * 256], 54);

    // 4. NaN FOV rejection, zero/negative binsize & min_fracexp safety, missing weights safety (#47-3)
    assert!(validate_fov(f64::NAN, 0.0).is_err());
    assert!(validate_fov(0.0, f64::INFINITY).is_err());
    let mut bcfg = BinDataConfig {
        time_binsize: 0.0,
        min_fracexp: -1.0,
        mask_weight: true,
        source_coords: Some((0.0, 0.0)),
        ..Default::default()
    };
    assert_eq!(bcfg.resolve_source_radec(83.63, 22.01), (0.0, 0.0));
    bcfg.source_coords = None;
    assert_eq!(bcfg.resolve_source_radec(83.63, 22.01), (83.63, 22.01));
    let mut ev1 = QuadEvents::new(0);
    ev1.time.push(100.0);
    ev1.cztseccnt.push(100.0);
    ev1.cztntick.push(0);
    ev1.pha.push(100);
    ev1.detid.push(0);
    ev1.pixid.push(0);
    ev1.detx.push(0);
    ev1.dety.push(0);
    ev1.veto.push(0);
    ev1.alpha.push(0);
    ev1.pi.push(100);
    ev1.energy.push(55.0);
    let (_s, lc, _w) = bindata_quad(
        &ev1,
        &BadpixTable::default(),
        &QuadLivetime::default(),
        None,
        &QuadGainCorCaldb::default(),
        true,
        &[0.0; 512],
        100.0,
        105.0,
        0.0,
        0.0,
        "NONE".into(),
        &bcfg,
    );
    assert!(lc.rate.iter().all(|v| v.is_finite()));

    let empty_rsp_caldb = QuadRspCaldb {
        n_ebins: 0,
        n_pibins: 0,
        emin: 1.0,
        emax: 301.0,
        pixresp_matrix: Vec::new(),
        resp_id: Vec::new(),
        eff_area: Vec::new(),
        eff_area_cor1: Vec::new(),
        eff_area_cor2: Vec::new(),
        lld: Vec::new(),
        uld: Vec::new(),
    };
    let empty_rsp = generate_quad_response(
        0,
        &empty_rsp_caldb,
        &BadpixTable::default(),
        0,
        &[],
        None,
        &[],
        Vec::new(),
        Vec::new(),
    );
    assert_eq!(empty_rsp.n_ebins, 0);

    // 5. eventsep drops 3+ event bursts from both singles & doubles; flickpix handles huge time gaps (#47-2)
    let mut ev_tri = QuadEvents::new(0);
    for (i, &t) in [1.0, 2.0, 2.000_010, 2.000_020, 3.0, 3.000_020, 1.0e11].iter().enumerate() {
        ev_tri.time.push(t);
        ev_tri.cztseccnt.push(t);
        ev_tri.cztntick.push(0);
        ev_tri.pha.push(100);
        ev_tri.detid.push(0);
        ev_tri.pixid.push(i as u8);
        ev_tri.detx.push(i as u8);
        ev_tri.dety.push(63);
        ev_tri.veto.push(0);
        ev_tri.alpha.push(0);
        ev_tri.pi.push(100);
        ev_tri.energy.push(55.0);
    }
    let (s_ev, d_ev) = eventsep_quad(&ev_tri, 30e-6);
    assert_eq!(s_ev.len(), 2, "3-event burst must be excluded from singles");
    assert_eq!(d_ev.len(), 2, "3-event burst must be excluded from doubles");
    let fp = flickpixclean_quad(
        &ev_tri,
        &vec![0u8; PIXELS_PER_QUAD],
        &vec![1.0; PIXELS_PER_QUAD],
        &NoiseCleanConfig::default(),
    );
    assert_eq!(fp.events.len(), ev_tri.len());
    assert!(pix_threshold_fix(100, 3500) >= 2);

    // 6. Undersized image arrays never panic & sky image WCS writes TX---TAN / TY---TAN (#47-4, #55-4)
    let tiny = vec![vec![1.0f32; 4]; 4];
    let asp2 = QuadAspect { quad_id: 2, ..Default::default() };
    let img_b = balanced_cross_correlation_from_mask(&tiny, &tiny, 2, &asp2);
    let img_f = fft_cross_correlation(&tiny, &tiny, 1, &asp2);
    assert_eq!(img_b.width, 5);
    assert!(img_f.data.iter().all(|v| v.is_finite() && *v >= 0.0));
    let dir = tempdir().unwrap();
    let img_path = dir.path().join("wcs_test.img");
    write_sky_image_file(&img_path, "Q2", &CztiObsMeta::default(), &img_f).unwrap();
    let img_fits = FitsFile::open(&img_path).unwrap();
    let h = img_fits.hdu_by_name("Q2").unwrap().header();
    assert_eq!(h.get_str("CTYPE1B").map(str::trim), Some("TX---TAN"));
    assert_eq!(h.get_str("CTYPE2B").map(str::trim), Some("TY---TAN"));

    // 7. run_addspec (`cztaddspec`) spectrum + exposure-weighted RMF co-addition (#55-7)
    let meta = CztiObsMeta { tstart: 100.0, tstop: 200.0, ..Default::default() };
    let mk_spec = |exp: f64, c: f32| CztiSpectrum {
        quad_id: Some(0),
        channel: (0..512).collect(),
        counts: vec![c; 512],
        stat_err: vec![c.sqrt(); 512],
        sys_err: vec![0.02; 512],
        quality: vec![0; 512],
        grouping: vec![1; 512],
        exposure: exp,
        thetax_deg: 0.0,
        thetay_deg: 0.0,
        mask_weighted: true,
        badpix_threshold: 0,
        gti_type: GtiType::Quad,
        respfile: "r.rsp".into(),
        tstart: 100.0,
        tstop: 200.0,
    };
    let eb = EboundsTable::default();
    let mk_rsp = |val: f32| CztiResponse {
        quad_id: Some(0),
        n_ebins: 2,
        n_pibins: 512,
        energ_lo: vec![10.0, 20.0],
        energ_hi: vec![20.0, 30.0],
        matrix: {
            let mut m = vec![0.0f32; 2 * 512];
            m[10] = val;
            m[512 + 20] = val;
            m
        },
        ebounds_channel: eb.channel.clone(),
        ebounds_emin: eb.e_min.clone(),
        ebounds_emax: eb.e_max.clone(),
    };
    let p1 = dir.path().join("s1.pha");
    let p2 = dir.path().join("s2.pha");
    let r1 = dir.path().join("r1.rsp");
    let r2 = dir.path().join("r2.rsp");
    let po = dir.path().join("out.pha");
    let ro = dir.path().join("out.rsp");
    write_spectrum_file(&p1, &meta, &mk_spec(100.0, 10.0)).unwrap();
    write_spectrum_file(&p2, &meta, &mk_spec(300.0, 30.0)).unwrap();
    write_response_file(&r1, &meta, &mk_rsp(100.0)).unwrap();
    write_response_file(&r2, &meta, &mk_rsp(200.0)).unwrap();
    let (cs, cr) = astrosat_czti::run_addspec(&[p1, p2], &[r1, r2], &po, Some(&ro), true).unwrap();
    assert!((cs.exposure - 400.0).abs() < 1e-6);
    assert!((cs.counts[0] - 40.0).abs() < 1e-4);
    // Exposure-weighted response: (100*100 + 300*200) / 400 = 175.0
    assert!((cr.unwrap().matrix[10] - 175.0).abs() < 1e-3);

    // 8. Neon review regressions (#101 + #103):
    // 8a. flickpixclean temporal per-pixel sigma (steady bright pixel NOT flagged; transient spike FLAGGED; >=120ks uses 8-sigma)
    let mut ev_fp = QuadEvents::new(0);
    ev_fp.exposure = 5_000.0; // <20ks -> (5.0 sigma, 1 bin)
    for bin in 0..40usize {
        let tb = 1_000.0 + (bin as f64) * 100.0;
        // 10 steady background pixels (1 event each) + 1 steady bright Crab pixel (pix 0: 10 events EVERY bin)
        for p in 1..=10u8 {
            ev_fp.time.push(tb + 1.0 + (p as f64) * 0.5);
            ev_fp.cztseccnt.push(tb + 1.0);
            ev_fp.cztntick.push(0);
            ev_fp.pha.push(50);
            ev_fp.detid.push(0);
            ev_fp.pixid.push(p);
            ev_fp.detx.push(p);
            ev_fp.dety.push(63);
            ev_fp.veto.push(0);
            ev_fp.alpha.push(0);
            ev_fp.pi.push(40);
            ev_fp.energy.push(25.0);
        }
        for k in 0..10 {
            ev_fp.time.push(tb + 10.0 + (k as f64) * 5.0);
            ev_fp.cztseccnt.push(tb + 10.0);
            ev_fp.cztntick.push(0);
            ev_fp.pha.push(50);
            ev_fp.detid.push(0);
            ev_fp.pixid.push(0);
            ev_fp.detx.push(0);
            ev_fp.dety.push(63);
            ev_fp.veto.push(0);
            ev_fp.alpha.push(0);
            ev_fp.pi.push(40);
            ev_fp.energy.push(25.0);
        }
        // Flickering pixel (pix 42): 0 events in 39 bins, 20 events in bin 20 (>6 sigma temporal spike!)
        if bin == 20 {
            for k in 0..20 {
                ev_fp.time.push(tb + 60.0 + (k as f64) * 1.5);
                ev_fp.cztseccnt.push(tb + 60.0);
                ev_fp.cztntick.push(0);
                ev_fp.pha.push(50);
                ev_fp.detid.push(0);
                ev_fp.pixid.push(42);
                ev_fp.detx.push(10);
                ev_fp.dety.push(61);
                ev_fp.veto.push(0);
                ev_fp.alpha.push(0);
                ev_fp.pi.push(40);
                ev_fp.energy.push(25.0);
            }
        }
    }
    let bp_init = vec![pix_flag::GOOD; PIXELS_PER_QUAD];
    let pe_init = vec![1.0f64; PIXELS_PER_QUAD];
    let gti_fp = vec![(1_000.0, 5_000.0)];
    let fp_res = flickpixclean_quad_with_gti(
        &ev_fp,
        &bp_init,
        &pe_init,
        Some(&gti_fp),
        &NoiseCleanConfig::default(),
    );
    assert_eq!(
        fp_res.badpix_flags[0],
        pix_flag::GOOD,
        "steadily illuminated bright pixel must NOT be falsely flagged as flickering"
    );
    assert_eq!(
        fp_res.badpix_flags[42],
        pix_flag::FLICKERING,
        "temporally flickering pixel must be flagged as FLICKERING"
    );
    // For >=120 ks exposure, threshold rises to 8-sigma / 2 bins, so single-bin spike is NOT flagged
    let mut ev_fp_long = ev_fp.clone();
    ev_fp_long.exposure = 150_000.0;
    let fp_long = flickpixclean_quad_with_gti(
        &ev_fp_long,
        &bp_init,
        &pe_init,
        Some(&gti_fp),
        &NoiseCleanConfig::default(),
    );
    assert_eq!(
        fp_long.badpix_flags[42],
        pix_flag::GOOD,
        "single spike bin must not trigger >=120ks (8-sigma, 2-bin) threshold"
    );

    // 8b. superbunchclean pairwise 3x3 neighborhood guard: adjacent mag==1 pairs (d=1.0) MUST count
    let mut ev_sb = QuadEvents::new(0);
    // Create a dense 2x4 cluster of 8 single-count pixels (mag=1 each, adjacent d=1) inside 0.1s
    // after a heavy bunch at t=100.0s: allowable_hot = sum(cij)/8 >> 3.0 -> triggers BTI!
    for dx in 0..4u8 {
        for dy in 0..2u8 {
            let t = 100.01 + ((dx * 2 + dy) as f64) * 0.005;
            ev_sb.time.push(t);
            ev_sb.cztseccnt.push(t);
            ev_sb.cztntick.push(0);
            ev_sb.pha.push(50);
            ev_sb.detid.push(0);
            ev_sb.pixid.push(dx + dy * 16);
            ev_sb.detx.push(dx);
            ev_sb.dety.push(dy);
            ev_sb.veto.push(0);
            ev_sb.alpha.push(0);
        }
    }
    let mut b_sb = BunchTable::default();
    b_sb.time.push(100.0);
    b_sb.time_dfs.push(0);
    b_sb.time_dsl.push(0);
    b_sb.num_event.push(30); // heavy bunch (>=15)
    b_sb.detid1.push(0);
    b_sb.detid_fevt.push(0);
    let lv_sb = QuadLivetime {
        quad_id: 0,
        time: (95..110).map(|s| s as f64 + 0.5).collect(),
        fracexp: vec![1.0; 15],
        binsize: 1.0,
    };
    let cfg_sb = NoiseCleanConfig {
        allowable_hot_thresh: 1.0,
        ..Default::default()
    };
    let sb_res = superbunchclean_quad(&ev_sb, &b_sb, &lv_sb, &cfg_sb);
    assert!(
        !sb_res.bti_intervals.is_empty(),
        "adjacent mag==1 cluster must satisfy pairwise counter4>2 and generate heavy-bunch BTI"
    );

    // 8c. bindata non-dyadic frac_binsize remainder wrap-around (`0.1 + 3e-8` with `time_binsize=1.0`)
    let mut ev_bd = QuadEvents::new(0);
    ev_bd.tstart = 200_000_000.0;
    ev_bd.tstop = 200_000_002.0;
    ev_bd.exposure = 1.0;
    for k in 0..10 {
        let t = 200_000_000.05 + (k as f64) * 0.1;
        ev_bd.time.push(t);
        ev_bd.cztseccnt.push(t);
        ev_bd.cztntick.push(0);
        ev_bd.pha.push(100);
        ev_bd.detid.push(0);
        ev_bd.pixid.push(0);
        ev_bd.detx.push(0);
        ev_bd.dety.push(63);
        ev_bd.veto.push(0);
        ev_bd.alpha.push(0);
        ev_bd.pi.push(100);
        ev_bd.energy.push(55.0);
    }
    let lv_nd = QuadLivetime {
        quad_id: 0,
        // 0.1 s bins with MET f64 ULP perturbation (~3e-8 s) where `1.0 % frac_binsize ≈ 0.0999997`
        time: (0..25)
            .map(|k| 200_000_000.05 + (k as f64) * (0.1 + 3e-8))
            .collect(),
        fracexp: vec![0.05; 25], // 50% livetime in each 0.1s bin -> 0.5 per 1.0s bin
        binsize: 0.1,
    };
    let bd_cfg = BinDataConfig {
        mask_weight: false,
        time_binsize: 1.0,
        ..Default::default()
    };
    let (_sp_nd, lc_nd, _) = bindata_quad(
        &ev_bd,
        &BadpixTable::default(),
        &lv_nd,
        None,
        &QuadGainCorCaldb::default(),
        false,
        &[0.0; PI_CHANNELS],
        200_000_000.0,
        200_000_002.0,
        0.0,
        0.0,
        "NONE".into(),
        &bd_cfg,
    );
    assert!(
        !lc_nd.fracexp.is_empty() && (lc_nd.fracexp[0] - 0.5).abs() < 0.05,
        "non-dyadic livetime wrap-around must apply FRACEXP ~0.5 (not fall back to 1.0), got {:?}",
        lc_nd.fracexp
    );

    // 8d. read_response_file rejects corrupt huge DETCHANS cleanly without OOM/panic (#103)
    let bad_rsp_path = dir.path().join("bad_detchans.rsp");
    let mut bad_rsp_hdu = astrosat_fits::BinTableBuilder::new("MATRIX")
        .add_f32("ENERG_LO", Some("keV"), vec![10.0, 20.0])
        .add_f32("ENERG_HI", Some("keV"), vec![20.0, 30.0])
        .build()
        .unwrap();
    bad_rsp_hdu.header.set_i64("DETCHANS", 500_000_000, "corrupt huge channel count");
    let mut bad_hdus = [astrosat_fits::Hdu::empty_primary(), bad_rsp_hdu];
    astrosat_fits::write_fits(&bad_rsp_path, &mut bad_hdus, true).unwrap();
    assert!(
        read_response_file(&bad_rsp_path).is_err(),
        "read_response_file must reject huge DETCHANS cleanly with Err"
    );
}
