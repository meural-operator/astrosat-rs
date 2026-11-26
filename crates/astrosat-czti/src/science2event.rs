//! Stage 1a: `cztscience2event` — Level-1 2048-byte packet decoding & event generation.
//!
//! Decodes AstroSat CZTI Level-1 `modeM0` and `modeSS` binary packets (2048 bytes / 1024 big-endian
//! 16-bit words) into calibrated Level-2 event lists (`Q0..Q3`), cosmic-ray onboard bunch tables,
//! `VETOSPECTRUM`, `SSM Data`, `TEMP`, and frame header housekeeping (`*.hdr`).

use std::collections::HashMap;
use serde::{Deserialize, Serialize};

use crate::types::{
    detid_pixid_to_detxy, BunchTable, FrameHeaderTable, QuadEvents, SsmDataTable, TempTable,
    VetoSpectrumTable, CZTSPEC_SIZE, CZT_TICK_HZ, NUM_DET_PER_QUAD, VETOSPEC_SIZE,
};

/// Size of a single CZTI Level-1 telemetry packet in bytes.
pub const PACKET_SIZE_BYTES: usize = 2048;
/// Number of 16-bit words in a single CZTI Level-1 telemetry packet.
pub const PACKET_WORDS: usize = 1024;
/// Packet header synchronization word 1.
pub const SYNC1: u16 = 0xF9A4;
/// Packet header synchronization word 2.
pub const SYNC2: u16 = 0x2BB1;

/// Configuration for `science2event`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Science2EventConfig {
    /// Input packet byte order (`true` for standard AstroSat big-endian Level-1 packets).
    pub big_endian: bool,
}

impl Default for Science2EventConfig {
    fn default() -> Self {
        Self { big_endian: true }
    }
}

/// Decoded Front-End Board (FEB) header from packet 0 of a CZTI frame.
#[derive(Debug, Clone, Default)]
pub struct FebHeader {
    pub veto_spec_range: bool,
    pub base_add: u8,
    pub buffer_no: bool,
    pub cmd_status: bool,
    pub event_read_mode: bool,
    pub czt_status: bool,
    pub czt_no: u8,
    pub channel_no: u8,
    pub adc_output: u16,
    pub cmd_1sec: u16,
    pub alpha_count: u16,
    pub veto_count: u16,
    pub czt_count_lt_uld: u16,
    pub czt_count_gt_uld: u16,
    pub czt_data_read: u16,
    pub temperature: [u16; NUM_DET_PER_QUAD],
    pub has_temperature: bool,
}

impl FebHeader {
    /// Decode the FEB header from packet 0 words (`data` length = 1024).
    pub fn decode(words: &[u16; PACKET_WORDS], mode_id: u8, data_id: u8) -> Option<Self> {
        let word_feb = if (4..=11).contains(&data_id) {
            24usize
        } else if mode_id <= 7 && data_id <= 3 {
            8usize
        } else {
            return None;
        };

        let s = 16usize;
        let w0 = words[s];
        let w1 = words[s + 1];
        let mut feb = Self {
            veto_spec_range: ((w0 >> 11) & 0x1) != 0,
            base_add: ((w0 >> 8) & 0x7) as u8,
            buffer_no: ((w0 >> 7) & 0x1) != 0,
            cmd_status: ((w0 >> 6) & 0x1) != 0,
            event_read_mode: ((w0 >> 5) & 0x1) != 0,
            czt_status: ((w0 >> 4) & 0x1) != 0,
            czt_no: (w0 & 0xF) as u8,
            channel_no: ((w1 >> 12) & 0x7) as u8,
            adc_output: w1 & 0x0FFF,
            cmd_1sec: words[s + 2],
            alpha_count: words[s + 3],
            veto_count: words[s + 4],
            czt_count_lt_uld: words[s + 5],
            czt_count_gt_uld: words[s + 6],
            czt_data_read: words[s + 7],
            temperature: [0u16; NUM_DET_PER_QUAD],
            has_temperature: word_feb == 24,
        };

        if word_feb == 24 {
            for i in 0..NUM_DET_PER_QUAD {
                feb.temperature[i] = words[s + 8 + i];
            }
        }
        Some(feb)
    }
}

/// Decoded single 2048-byte telemetry packet.
#[derive(Debug, Clone)]
pub struct DecodedPacket {
    pub sync1: u16,
    pub sync2: u16,
    pub mml: u8,
    pub data_id: u8,
    pub pkt_no: u8,
    pub mode_id: u8,
    pub quad_id: i8,
    pub wcnt: u16,
    pub is_first: bool,
    pub is_valid: bool,
    pub frame_header: Option<PacketFrameHeader>,
    pub feb_header: Option<FebHeader>,
    pub hbt_data: Vec<u16>,
}

/// Frame-level header fields present only in packet 0 (`pkt_no == 0`).
#[derive(Debug, Clone, Default)]
pub struct PacketFrameHeader {
    pub frame_no: u16,
    pub status_flag: u16,
    pub write_pkt_no: u16,
    pub read_pkt_no: u16,
    pub command: i32,
    pub second_count: u32,
    pub laxpc_time: f64,
    pub dcnt: u16,
    pub boot_page_no: u8,
    pub error_count: u8,
    pub error_flag: u8,
}

/// Assembled CZTI telemetry frame consisting of packet 0 followed by continuation packets.
#[derive(Debug, Clone)]
pub struct DecodedFrame {
    pub data_id: u8,
    pub mode_id: u8,
    pub quad_id: i8,
    pub header: PacketFrameHeader,
    pub feb: FebHeader,
    pub packets: Vec<DecodedPacket>,
}

/// Result of decoding one quadrant's packet stream (`CZT_QUAD1..4`).
#[derive(Debug, Clone, Default)]
pub struct DecodedQuadStream {
    pub events: QuadEvents,
    pub bunches: BunchTable,
    pub veto_spec: VetoSpectrumTable,
    pub ssm_data: SsmDataTable,
    pub temp: TempTable,
    pub headers: FrameHeaderTable,
    pub total_packets: usize,
    pub valid_packets: usize,
    pub invalid_packets: usize,
    pub total_frames: usize,
    pub missing_frames: usize,
}

/// Convert a 2048-byte raw packet slice into 1024 16-bit words.
#[inline]
pub fn packet_bytes_to_words(pkt: &[u8], big_endian: bool) -> [u16; PACKET_WORDS] {
    debug_assert!(pkt.len() >= PACKET_SIZE_BYTES);
    let mut words = [0u16; PACKET_WORDS];
    if big_endian {
        for (dst, &chunk) in words.iter_mut().zip(pkt.as_chunks::<2>().0) {
            *dst = u16::from_be_bytes(chunk);
        }
    } else {
        for (dst, &chunk) in words.iter_mut().zip(pkt.as_chunks::<2>().0) {
            *dst = u16::from_le_bytes(chunk);
        }
    }
    words
}

/// Decode a single 2048-byte CZTI packet following `l1evtdecode.cpp::L1evtData::decode`.
pub fn decode_packet(pkt: &[u8], big_endian: bool) -> DecodedPacket {
    if pkt.len() < PACKET_SIZE_BYTES {
        return DecodedPacket {
            sync1: 0,
            sync2: 0,
            mml: 0,
            data_id: 0,
            pkt_no: 0,
            mode_id: 0,
            quad_id: -1,
            wcnt: 0,
            is_first: false,
            is_valid: false,
            frame_header: None,
            feb_header: None,
            hbt_data: Vec::new(),
        };
    }

    let w = packet_bytes_to_words(pkt, big_endian);
    let sync1 = w[0];
    let sync2 = w[1];
    let mml = ((w[2] >> 13) & 0x3) as u8;
    let data_id = ((w[2] >> 8) & 0x1F) as u8;
    let pkt_no = ((w[2] >> 4) & 0xF) as u8;
    let mode_id = (w[2] & 0xF) as u8;
    let raw_wcnt = w[3];
    let is_first = pkt_no == 0;

    let (wcnt, frame_header) = if is_first {
        let command = (((w[9] as u32) << 16) | (w[8] as u32)) as i32;
        let second_count = ((w[10] as u32) << 16) | (w[11] as u32);
        let laxpc_sec = (((w[14] as u32) << 12) & 0xFF000) | ((w[13] as u32) >> 4);
        let laxpc_ms = ((w[12] >> 10) | ((w[13] & 0xF) << 6)) as f64;
        let laxpc_us = (w[12] & 0x03FF) as f64;
        let laxpc_time = (laxpc_sec as f64) + laxpc_ms / 1000.0 + laxpc_us / 1_000_000.0;
        let fh = PacketFrameHeader {
            frame_no: w[4],
            status_flag: w[5],
            write_pkt_no: w[6],
            read_pkt_no: w[7],
            command,
            second_count,
            laxpc_time,
            dcnt: w[15] & 0x3FFF,
            boot_page_no: ((w[15] >> 14) & 0x3) as u8,
            error_count: ((w[14] >> 12) & 0xF) as u8,
            error_flag: ((w[14] >> 8) & 0xF) as u8,
        };
        let adjusted = raw_wcnt.checked_sub(12);
        (adjusted, Some(fh))
    } else {
        (Some(raw_wcnt), None)
    };

    let quad_id: i8 = if data_id < 12 { (data_id % 4) as i8 } else { -1 };

    let valid_wcnt = match wcnt {
        Some(wc) if is_first => (wc as usize) <= (PACKET_WORDS - 16),
        Some(wc) => (wc as usize) <= (PACKET_WORDS - 4),
        None => false,
    };

    let is_valid = sync1 == SYNC1
        && sync2 == SYNC2
        && data_id < 16
        && mml < 4
        && mode_id < 16
        && valid_wcnt;

    let (feb_header, hbt_data) = if is_valid {
        let wc = wcnt.unwrap_or(0) as usize;
        let start = if is_first { 16 } else { 4 };
        let feb = if is_first {
            FebHeader::decode(&w, mode_id, data_id)
        } else {
            None
        };
        let data = w[start..start + wc].to_vec();
        (feb, data)
    } else {
        (None, Vec::new())
    };

    DecodedPacket {
        sync1,
        sync2,
        mml,
        data_id,
        pkt_no,
        mode_id,
        quad_id,
        wcnt: wcnt.unwrap_or(0),
        is_first,
        is_valid,
        frame_header,
        feb_header,
        hbt_data,
    }
}

/// Assemble raw 2048-byte packets into contiguous valid frames (mirroring `segregate_data`).
pub fn assemble_frames(raw_packets: &[&[u8]], big_endian: bool) -> (Vec<DecodedFrame>, usize, usize) {
    let mut frames: Vec<DecodedFrame> = Vec::new();
    let mut current_valid = false;
    let mut invalid_packets = 0usize;

    for raw in raw_packets {
        let pkt = decode_packet(raw, big_endian);
        if !pkt.is_valid {
            invalid_packets += 1;
            continue;
        }
        if pkt.is_first {
            current_valid = true;
            let header = pkt.frame_header.clone().unwrap_or_default();
            let feb = pkt.feb_header.clone().unwrap_or_default();
            let data_id = pkt.data_id;
            let mode_id = pkt.mode_id;
            let quad_id = pkt.quad_id;
            frames.push(DecodedFrame {
                data_id,
                mode_id,
                quad_id,
                header,
                feb,
                packets: vec![pkt],
            });
        } else if let Some(cur) = frames.last_mut() {
            if let Some(last_pkt) = cur.packets.last() {
                if pkt.pkt_no < last_pkt.pkt_no {
                    current_valid = false;
                }
            }
            if current_valid {
                cur.packets.push(pkt);
            } else {
                invalid_packets += 1;
            }
        } else {
            // Non-zero packet number before first pktNo==0
            invalid_packets += 1;
        }
    }

    let valid_packets: usize = frames.iter().map(|f| f.packets.len()).sum();
    (frames, valid_packets, invalid_packets)
}

/// Remove zero or junk (`> 50 C`) module temperatures in-place, matching
/// `TemperatureData::remove_junk_temperature` in `l1evtdecode.cpp`.
pub fn clean_junk_temperatures(temps: &mut [[f32; NUM_DET_PER_QUAD]]) {
    for row in temps.iter_mut() {
        if row[0] == 0.0 || row[0] > 50.0 {
            // Compute mode of valid detector temperatures in this row
            let mut counts: HashMap<i32, usize> = HashMap::new();
            for &t in row.iter() {
                if t > 0.0 && t <= 50.0 {
                    *counts.entry(t.round() as i32).or_insert(0) += 1;
                }
            }
            let mode_t = counts
                .into_iter()
                .max_by_key(|&(_, c)| c)
                .map(|(k, _)| k as f32)
                .unwrap_or(10.0);
            row[0] = mode_t;
        }
        for idet in 1..NUM_DET_PER_QUAD {
            if row[idet] == 0.0 || row[idet] > 50.0 {
                row[idet] = row[idet - 1];
            }
        }
    }
}

/// Buffered 3-record onboard bunch accumulator so incomplete bunches at a frame boundary
/// never push partial rows into `QuadEvents` or `BunchTable`.
struct PendingOnboardBunch {
    time_dfs: u8,
    time_dsl: u8,
    num_evts: u8,
    // Record 0 (first event)
    r0_pha: u16,
    r0_detid: u8,
    r0_pixid: u8,
    r0_veto: u16,
    r0_alpha: u8,
    r0_detx: u8,
    r0_dety: u8,
    // Record 1 (second event + timestamps for all 3 records)
    t1_czt: f64,
    t1_ut: f64,
    t1_ntick: u16,
    t2_czt: f64,
    t2_ut: f64,
    t2_ntick: u16,
    t3_czt: f64,
    t3_ut: f64,
    t3_ntick: u16,
    r1_pha: u16,
    r1_detid: u8,
    r1_pixid: u8,
    r1_veto: u16,
    r1_alpha: u8,
    r1_detx: u8,
    r1_dety: u8,
}

impl Default for PendingOnboardBunch {
    fn default() -> Self {
        Self {
            time_dfs: 0,
            time_dsl: 0,
            num_evts: 0,
            r0_pha: 0,
            r0_detid: 0,
            r0_pixid: 0,
            r0_veto: 0,
            r0_alpha: 0,
            r0_detx: 0,
            r0_dety: 0,
            t1_czt: 0.0,
            t1_ut: 0.0,
            t1_ntick: 0,
            t2_czt: 0.0,
            t2_ut: 0.0,
            t2_ntick: 0,
            t3_czt: 0.0,
            t3_ut: 0.0,
            t3_ntick: 0,
            r1_pha: 0,
            r1_detid: 0,
            r1_pixid: 0,
            r1_veto: 0,
            r1_alpha: 0,
            r1_detx: 0,
            r1_dety: 0,
        }
    }
}

/// Decode a quadrant's raw Level-1 packet buffer (`nrows * 2048` bytes) using a TCT time converter.
///
/// `tct_convert` maps instrument `cztseccnt` (seconds) to AstroSat UTC `Time` (seconds since MJD 55197).
pub fn decode_quad_packets<F>(
    quad_id: u8,
    raw_bytes: &[u8],
    cfg: &Science2EventConfig,
    tct_convert: F,
) -> DecodedQuadStream
where
    F: Fn(f64) -> f64,
{
    let n_pkts = raw_bytes.len() / PACKET_SIZE_BYTES;
    let pkt_slices: Vec<&[u8]> = (0..n_pkts)
        .map(|i| &raw_bytes[i * PACKET_SIZE_BYTES..(i + 1) * PACKET_SIZE_BYTES])
        .collect();

    let (frames, valid_packets, invalid_packets) = assemble_frames(&pkt_slices, cfg.big_endian);

    let mut missing_frames = 0usize;
    for pair in frames.windows(2) {
        let diff = (pair[1].header.frame_no as i32) - (pair[0].header.frame_no as i32);
        if diff > 1 {
            missing_frames += (diff - 1) as usize;
        }
    }

    let mut out = DecodedQuadStream {
        events: QuadEvents::new(quad_id),
        bunches: BunchTable::new(quad_id),
        veto_spec: VetoSpectrumTable::default(),
        ssm_data: SsmDataTable::default(),
        temp: TempTable::default(),
        headers: FrameHeaderTable::default(),
        total_packets: n_pkts,
        valid_packets,
        invalid_packets,
        total_frames: frames.len(),
        missing_frames,
    };

    for frame in &frames {
        let (feb_veto_size, eventsize, vetospec_start): (usize, usize, Option<usize>) =
            match frame.mode_id {
                0 | 4 => (240, 3, Some(8)),
                1 | 5 => (8, 3, None),
                2 | 6 => (244, 2, Some(8)),
                3 | 7 => (8, 2, None),
                _ => (240, 3, Some(8)),
            };

        let frame_sec = frame.header.second_count as f64;
        let frame_time = tct_convert(frame_sec);

        // Record frame header row
        let tot_words: usize = frame.packets.iter().map(|p| p.wcnt as usize).sum();
        let evt_word_cnt = tot_words.saturating_sub(feb_veto_size) / eventsize;
        out.headers.time.push(frame_time);
        out.headers.cztseccnt.push(frame_sec);
        out.headers.data_id.push(frame.data_id);
        out.headers.mode_id.push(frame.mode_id);
        out.headers.word_count.push(evt_word_cnt as i32);
        out.headers.frame_no.push(frame.header.frame_no);
        out.headers.status_flag.push(frame.header.status_flag);
        out.headers.write_pkt_no.push(frame.header.write_pkt_no);
        out.headers.read_pkt_no.push(frame.header.read_pkt_no);
        out.headers.command.push(frame.header.command);
        out.headers.laxpc_time.push(frame.header.laxpc_time);
        out.headers.dcnt.push(frame.header.dcnt);
        out.headers.error_count.push(frame.header.error_count);
        out.headers.error_flag.push(frame.header.error_flag);
        out.headers.boot_page_no.push(frame.header.boot_page_no);
        out.headers.czt_no.push(frame.feb.czt_no);
        out.headers.czt_status.push(frame.feb.czt_status as u8);
        out.headers.evt_read_mode.push(frame.feb.event_read_mode as u8);
        out.headers.cmd_status.push(frame.feb.cmd_status as u8);
        out.headers.buffer_no.push(frame.feb.buffer_no as u8);
        out.headers.base_add.push(frame.feb.base_add);
        out.headers.veto_spec_range.push(frame.feb.veto_spec_range as u8);
        out.headers.channel_no.push(frame.feb.channel_no);
        out.headers.adc_output.push(frame.feb.adc_output);
        out.headers.cmd_1sec.push(frame.feb.cmd_1sec);
        out.headers.veto_count.push(frame.feb.veto_count);
        out.headers.alpha_count.push(frame.feb.alpha_count);
        out.headers.czt_count_lt_uld.push(frame.feb.czt_count_lt_uld);
        out.headers.czt_count_gt_uld.push(frame.feb.czt_count_gt_uld);
        out.headers.czt_data_read.push(frame.feb.czt_data_read);

        // Skip non-event frames for event decoding
        if frame.mode_id > 7
            || (4..=7).contains(&frame.data_id)
            || (12..=15).contains(&frame.data_id)
            || !frame.feb.event_read_mode
        {
            continue;
        }

        // Decode event records for dataID in 0..=3
        if frame.data_id <= 3 {
            let mut nframe_bunch: u8 = 0;
            let mut bunch_sync: u8 = 0;
            let mut pending = PendingOnboardBunch::default();

            for pkt in &frame.packets {
                let (start_idx, n_events) = if pkt.pkt_no == 0 {
                    if pkt.hbt_data.len() <= feb_veto_size {
                        (0, 0)
                    } else {
                        (feb_veto_size, (pkt.hbt_data.len() - feb_veto_size) / eventsize)
                    }
                } else {
                    (0, pkt.hbt_data.len() / eventsize)
                };

                let mut i_evt = 0usize;
                let mut idx = start_idx;
                while i_evt < n_events {
                    if eventsize == 3 {
                        if nframe_bunch == 0 {
                            bunch_sync = ((pkt.hbt_data[idx] >> 13) & 0x7) as u8;
                        }

                        if bunch_sync != 7 {
                            let cztntick = pkt.hbt_data[idx];
                            let cztseccnt = (cztntick as f64) / CZT_TICK_HZ + frame_sec;
                            let pha = ((pkt.hbt_data[idx + 1] >> 4) & 0x0FFF) >> 2;
                            let detid = (pkt.hbt_data[idx + 1] & 0xF) as u8;
                            let pixid = ((pkt.hbt_data[idx + 2] >> 8) & 0xFF) as u8;
                            let veto = (pkt.hbt_data[idx + 2] >> 1) & 0x7F;
                            let alpha = (pkt.hbt_data[idx + 2] & 0x1) as u8;
                            let q = if pkt.quad_id >= 0 { pkt.quad_id as u8 } else { quad_id };
                            let (detx, dety) = detid_pixid_to_detxy(q, detid, pixid);

                            out.events.time.push(tct_convert(cztseccnt));
                            out.events.cztseccnt.push(cztseccnt);
                            out.events.cztntick.push(cztntick);
                            out.events.pha.push(pha);
                            out.events.detid.push(detid);
                            out.events.pixid.push(pixid);
                            out.events.detx.push(detx);
                            out.events.dety.push(dety);
                            out.events.veto.push(veto);
                            out.events.alpha.push(alpha);
                        } else if frame.header.boot_page_no == 2 {
                            // On-board bunch report spanning 3 consecutive 3-word records.
                            // Buffered locally in `pending` and committed atomically on record 2
                            // so a frame ending mid-bunch can never desynchronize event columns.
                            let q = if pkt.quad_id >= 0 { pkt.quad_id as u8 } else { quad_id };
                            if nframe_bunch == 0 {
                                pending.time_dfs = (pkt.hbt_data[idx] & 0x0001) as u8;
                                pending.time_dsl = ((pkt.hbt_data[idx] >> 1) & 0x3F) as u8;
                                pending.num_evts = ((pkt.hbt_data[idx] >> 7) & 0x3F) as u8;
                                pending.r0_pha = ((pkt.hbt_data[idx + 1] >> 4) & 0x0FFF) >> 2;
                                pending.r0_detid = (pkt.hbt_data[idx + 1] & 0xF) as u8;
                                pending.r0_pixid = ((pkt.hbt_data[idx + 2] >> 8) & 0xFF) as u8;
                                pending.r0_veto = (pkt.hbt_data[idx + 2] >> 1) & 0x7F;
                                pending.r0_alpha = (pkt.hbt_data[idx + 2] & 0x1) as u8;
                                let (detx, dety) =
                                    detid_pixid_to_detxy(q, pending.r0_detid, pending.r0_pixid);
                                pending.r0_detx = detx;
                                pending.r0_dety = dety;

                                nframe_bunch = 1;
                                i_evt += 1;
                                idx += eventsize;
                                if i_evt >= n_events {
                                    continue;
                                }
                            }
                            if nframe_bunch == 1 {
                                let cztntick = pkt.hbt_data[idx];
                                let cztseccnt = (cztntick as f64) / CZT_TICK_HZ + frame_sec;
                                let dfs = pending.time_dfs as u16;
                                let dsl = pending.time_dsl as u16;

                                pending.t1_czt = cztseccnt - (dfs as f64) / CZT_TICK_HZ;
                                pending.t1_ut = tct_convert(pending.t1_czt);
                                pending.t1_ntick = cztntick.saturating_sub(dfs);

                                pending.t2_czt = cztseccnt;
                                pending.t2_ut = tct_convert(cztseccnt);
                                pending.t2_ntick = cztntick;

                                pending.t3_czt = cztseccnt + (dsl as f64) / CZT_TICK_HZ;
                                pending.t3_ut = tct_convert(pending.t3_czt);
                                pending.t3_ntick = cztntick.saturating_add(dsl);

                                pending.r1_pha = ((pkt.hbt_data[idx + 1] >> 4) & 0x0FFF) >> 2;
                                pending.r1_detid = (pkt.hbt_data[idx + 1] & 0xF) as u8;
                                pending.r1_pixid = ((pkt.hbt_data[idx + 2] >> 8) & 0xFF) as u8;
                                pending.r1_veto = (pkt.hbt_data[idx + 2] >> 1) & 0x7F;
                                pending.r1_alpha = (pkt.hbt_data[idx + 2] & 0x1) as u8;
                                let (detx, dety) =
                                    detid_pixid_to_detxy(q, pending.r1_detid, pending.r1_pixid);
                                pending.r1_detx = detx;
                                pending.r1_dety = dety;

                                nframe_bunch = 2;
                                i_evt += 1;
                                idx += eventsize;
                                if i_evt >= n_events {
                                    continue;
                                }
                            }
                            if nframe_bunch == 2 {
                                let row_num = out.events.time.len() as i32;
                                let detid1 = (pkt.hbt_data[idx] & 0xF) as u8;
                                let detid2 = ((pkt.hbt_data[idx] >> 4) & 0xF) as u8;
                                let detid3 = ((pkt.hbt_data[idx] >> 8) & 0xF) as u8;
                                let detid4 = ((pkt.hbt_data[idx] >> 12) & 0xF) as u8;

                                let r2_pha = ((pkt.hbt_data[idx + 1] >> 4) & 0x0FFF) >> 2;
                                let r2_detid = (pkt.hbt_data[idx + 1] & 0xF) as u8;
                                let r2_pixid = ((pkt.hbt_data[idx + 2] >> 8) & 0xFF) as u8;
                                let r2_veto = (pkt.hbt_data[idx + 2] >> 1) & 0x7F;
                                let r2_alpha = (pkt.hbt_data[idx + 2] & 0x1) as u8;
                                let (r2_detx, r2_dety) = detid_pixid_to_detxy(q, r2_detid, r2_pixid);

                                // Atomically commit all 3 bunch events + 1 bunch table row
                                out.bunches.time.push(pending.t2_ut);
                                out.bunches.time_dfs.push(pending.time_dfs);
                                out.bunches.time_dsl.push(pending.time_dsl);
                                out.bunches.num_event.push(pending.num_evts);
                                out.bunches.detid1.push(detid1);
                                out.bunches.detid2.push(detid2);
                                out.bunches.detid3.push(detid3);
                                out.bunches.detid4.push(detid4);
                                out.bunches.row_num.push(row_num);
                                out.bunches.detid_fevt.push(pending.r0_detid);
                                out.bunches.pixid_fevt.push(pending.r0_pixid);
                                out.bunches.detid_sevt.push(pending.r1_detid);
                                out.bunches.pixid_sevt.push(pending.r1_pixid);
                                out.bunches.detid_levt.push(r2_detid);
                                out.bunches.pixid_levt.push(r2_pixid);

                                // Event 1 (first in bunch)
                                out.events.time.push(pending.t1_ut);
                                out.events.cztseccnt.push(pending.t1_czt);
                                out.events.cztntick.push(pending.t1_ntick);
                                out.events.pha.push(pending.r0_pha);
                                out.events.detid.push(pending.r0_detid);
                                out.events.pixid.push(pending.r0_pixid);
                                out.events.detx.push(pending.r0_detx);
                                out.events.dety.push(pending.r0_dety);
                                out.events.veto.push(pending.r0_veto);
                                out.events.alpha.push(pending.r0_alpha);

                                // Event 2 (second in bunch)
                                out.events.time.push(pending.t2_ut);
                                out.events.cztseccnt.push(pending.t2_czt);
                                out.events.cztntick.push(pending.t2_ntick);
                                out.events.pha.push(pending.r1_pha);
                                out.events.detid.push(pending.r1_detid);
                                out.events.pixid.push(pending.r1_pixid);
                                out.events.detx.push(pending.r1_detx);
                                out.events.dety.push(pending.r1_dety);
                                out.events.veto.push(pending.r1_veto);
                                out.events.alpha.push(pending.r1_alpha);

                                // Event 3 (last in bunch)
                                out.events.time.push(pending.t3_ut);
                                out.events.cztseccnt.push(pending.t3_czt);
                                out.events.cztntick.push(pending.t3_ntick);
                                out.events.pha.push(r2_pha);
                                out.events.detid.push(r2_detid);
                                out.events.pixid.push(r2_pixid);
                                out.events.detx.push(r2_detx);
                                out.events.dety.push(r2_dety);
                                out.events.veto.push(r2_veto);
                                out.events.alpha.push(r2_alpha);

                                nframe_bunch = 0;
                            }
                        }
                        // else: bunch_sync == 7 && boot_page_no != 2 -> corrupted event, ignored
                    } else if eventsize == 2 {
                        // Reduced mode 2-word event report (`modeID in {2,3,6,7}`)
                        let cztntick = (pkt.hbt_data[idx] & 0xFF80) >> 7;
                        let cztseccnt = (cztntick as f64) / CZT_TICK_HZ + frame_sec;
                        let pha = (((pkt.hbt_data[idx] & 0x7F) << 2)
                            | ((pkt.hbt_data[idx + 1] >> 14) & 0x3))
                            << 1;
                        let detid = (pkt.hbt_data[idx + 1] & 0xF) as u8;
                        let pixid = ((pkt.hbt_data[idx + 1] >> 6) & 0xFF) as u8;
                        let veto = (pkt.hbt_data[idx + 1] >> 5) & 0x1;
                        let alpha = ((pkt.hbt_data[idx + 1] >> 4) & 0x1) as u8;
                        let q = if pkt.quad_id >= 0 { pkt.quad_id as u8 } else { quad_id };
                        let (detx, dety) = detid_pixid_to_detxy(q, detid, pixid);

                        out.events.time.push(tct_convert(cztseccnt));
                        out.events.cztseccnt.push(cztseccnt);
                        out.events.cztntick.push(cztntick);
                        out.events.pha.push(pha);
                        out.events.detid.push(detid);
                        out.events.pixid.push(pixid);
                        out.events.detx.push(detx);
                        out.events.dety.push(dety);
                        out.events.veto.push(veto);
                        out.events.alpha.push(alpha);
                    }

                    i_evt += 1;
                    idx += eventsize;
                }
            }
        }

        // Extract Veto Spectrum if present (`vetospec_start == Some(8)`)
        if let Some(vstart) = vetospec_start {
            let mut vspec = [0u16; VETOSPEC_SIZE];
            let p0 = &frame.packets[0];
            if frame.feb.veto_spec_range {
                let n = VETOSPEC_SIZE.min(p0.hbt_data.len().saturating_sub(vstart));
                vspec[..n].copy_from_slice(&p0.hbt_data[vstart..vstart + n]);
            } else {
                let n = (VETOSPEC_SIZE - 24).min(p0.hbt_data.len().saturating_sub(vstart));
                vspec[24..24 + n].copy_from_slice(&p0.hbt_data[vstart..vstart + n]);
            }
            out.veto_spec.veto_spec.push(vspec);
            out.veto_spec.cztseccnt.push(frame_sec);
            out.veto_spec.time.push(frame_time);
            out.veto_spec
                .quad_id
                .push(if frame.quad_id >= 0 { frame.quad_id as u8 } else { quad_id });
        }

        // Extract SSM Data & Temperature (`dataID in 8..=11`)
        if (8..=11).contains(&frame.data_id) {
            let q = if frame.quad_id >= 0 { frame.quad_id as u8 } else { quad_id };
            if frame.packets.len() == 2 {
                let p0 = &frame.packets[0].hbt_data;
                let p1 = &frame.packets[1].hbt_data;
                if p0.len() >= 1008 && p1.len() >= 456 + CZTSPEC_SIZE {
                    let mut vs = [0u16; VETOSPEC_SIZE];
                    vs.copy_from_slice(&p0[184..184 + VETOSPEC_SIZE]);
                    let mut cs = [0u16; CZTSPEC_SIZE];
                    cs.copy_from_slice(&p0[440..440 + CZTSPEC_SIZE]);
                    let mut csv = [0u16; CZTSPEC_SIZE];
                    csv[..56].copy_from_slice(&p0[952..952 + 56]);
                    csv[56..].copy_from_slice(&p1[..456]);
                    let mut csa = [0u16; CZTSPEC_SIZE];
                    csa.copy_from_slice(&p1[456..456 + CZTSPEC_SIZE]);

                    out.ssm_data.time.push(frame_time);
                    out.ssm_data.cztseccnt.push(frame_sec);
                    out.ssm_data.veto_spec.push(vs);
                    out.ssm_data.czt_spec.push(cs);
                    out.ssm_data.czt_spec_veto.push(csv);
                    out.ssm_data.czt_spec_alpha.push(csa);
                    out.ssm_data.quad_id.push(q);
                }
            }

            // Temperature data (16 detector modules)
            let mut trow = [0.0f32; NUM_DET_PER_QUAD];
            for i in 0..NUM_DET_PER_QUAD {
                trow[i] = frame.feb.temperature[i] as f32;
            }
            out.temp.time.push(frame_time);
            out.temp.cztseccnt.push(frame_sec);
            out.temp.temperature.push(trow);
            out.temp.quad_flag.push(q);
        }
    }

    // Ensure truncated partial onboard bunch at frame boundary doesn't leave mismatched column lengths
    let n_evts = out.events.time.len();
    out.events.cztseccnt.truncate(n_evts);
    out.events.cztntick.truncate(n_evts);
    out.events.pha.truncate(n_evts);
    out.events.detid.truncate(n_evts);
    out.events.pixid.truncate(n_evts);
    out.events.detx.truncate(n_evts);
    out.events.dety.truncate(n_evts);
    out.events.veto.truncate(n_evts);
    out.events.alpha.truncate(n_evts);

    let n_bun = out
        .bunches
        .time
        .len()
        .min(out.bunches.detid4.len())
        .min(out.bunches.row_num.len());
    out.bunches.time.truncate(n_bun);
    out.bunches.time_dfs.truncate(n_bun);
    out.bunches.time_dsl.truncate(n_bun);
    out.bunches.num_event.truncate(n_bun);
    out.bunches.detid1.truncate(n_bun);
    out.bunches.detid2.truncate(n_bun);
    out.bunches.detid3.truncate(n_bun);
    out.bunches.detid4.truncate(n_bun);
    out.bunches.row_num.truncate(n_bun);
    out.bunches.detid_fevt.truncate(n_bun);
    out.bunches.pixid_fevt.truncate(n_bun);
    out.bunches.detid_sevt.truncate(n_bun);
    out.bunches.pixid_sevt.truncate(n_bun);
    out.bunches.detid_levt.truncate(n_bun);
    out.bunches.pixid_levt.truncate(n_bun);

    if !out.events.is_empty() {
        out.events.tstart = out.events.time[0];
        out.events.tstop = *out.events.time.last().unwrap_or(&0.0);
    }

    clean_junk_temperatures(&mut out.temp.temperature);

    out
}
