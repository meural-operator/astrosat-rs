//! FITS 4.0 header cards, values, `CONTINUE` long strings, and 2880-byte block serialization.

use crate::checksum::FITS_BLOCK_SIZE;
use crate::error::{FitsError, Result};

/// Value stored in a FITS header card.
#[derive(Debug, Clone, PartialEq)]
pub enum CardValue {
    Logical(bool),
    Integer(i64),
    Float(f64),
    Str(String),
    Empty,
}

/// A single logical FITS header record (may span multiple physical 80-byte cards via `CONTINUE`).
#[derive(Debug, Clone, PartialEq)]
pub struct Card {
    pub keyword: String,
    pub value: CardValue,
    pub comment: Option<String>,
}

impl Card {
    pub fn new_str(
        keyword: impl AsRef<str>,
        value: impl Into<String>,
        comment: impl Into<Option<String>>,
    ) -> Self {
        Self {
            keyword: normalize_keyword(keyword.as_ref()),
            value: CardValue::Str(value.into()),
            comment: comment.into(),
        }
    }

    pub fn new_int(
        keyword: impl AsRef<str>,
        value: i64,
        comment: impl Into<Option<String>>,
    ) -> Self {
        Self {
            keyword: normalize_keyword(keyword.as_ref()),
            value: CardValue::Integer(value),
            comment: comment.into(),
        }
    }

    pub fn new_float(
        keyword: impl AsRef<str>,
        value: f64,
        comment: impl Into<Option<String>>,
    ) -> Self {
        Self {
            keyword: normalize_keyword(keyword.as_ref()),
            value: CardValue::Float(value),
            comment: comment.into(),
        }
    }

    pub fn new_bool(
        keyword: impl AsRef<str>,
        value: bool,
        comment: impl Into<Option<String>>,
    ) -> Self {
        Self {
            keyword: normalize_keyword(keyword.as_ref()),
            value: CardValue::Logical(value),
            comment: comment.into(),
        }
    }

    pub fn comment(text: impl Into<String>) -> Self {
        Self {
            keyword: "COMMENT".to_string(),
            value: CardValue::Empty,
            comment: Some(text.into()),
        }
    }

    pub fn history(text: impl Into<String>) -> Self {
        Self {
            keyword: "HISTORY".to_string(),
            value: CardValue::Empty,
            comment: Some(text.into()),
        }
    }

    /// Format this logical card into one or more 80-byte physical FITS cards (handling `CONTINUE`).
    pub fn to_fits_cards(&self) -> Vec<[u8; 80]> {
        let kw = normalize_keyword(&self.keyword);

        // COMMENT / HISTORY / blank cards have no "= " indicator.
        if kw == "COMMENT" || kw == "HISTORY" || kw.is_empty() {
            let text = sanitize_ascii(self.comment.as_deref().unwrap_or(""));
            if text.is_empty() {
                let mut card = [b' '; 80];
                card[..kw.len().min(8)].copy_from_slice(&kw.as_bytes()[..kw.len().min(8)]);
                return vec![card];
            }
            let mut out = Vec::new();
            let bytes = text.as_bytes();
            let mut pos = 0;
            while pos < bytes.len() {
                let end = (pos + 72).min(bytes.len());
                let mut card = [b' '; 80];
                card[..kw.len().min(8)].copy_from_slice(&kw.as_bytes()[..kw.len().min(8)]);
                card[8..8 + (end - pos)].copy_from_slice(&bytes[pos..end]);
                out.push(card);
                pos = end;
            }
            return out;
        }

        match &self.value {
            CardValue::Empty => {
                let mut card = [b' '; 80];
                card[..kw.len().min(8)].copy_from_slice(&kw.as_bytes()[..kw.len().min(8)]);
                if let Some(c) = &self.comment {
                    let s = sanitize_ascii(c);
                    let slice = &s.as_bytes()[..s.len().min(72)];
                    card[8..8 + slice.len()].copy_from_slice(slice);
                }
                vec![card]
            }
            CardValue::Logical(b) => {
                let val_str = format!("{:>20}", if *b { "T" } else { "F" });
                vec![format_simple_card(&kw, &val_str, self.comment.as_deref())]
            }
            CardValue::Integer(n) => {
                let val_str = format!("{:>20}", n);
                vec![format_simple_card(&kw, &val_str, self.comment.as_deref())]
            }
            CardValue::Float(f) => {
                let formatted = format_fits_float(*f);
                let val_str = format!("{:>20}", formatted);
                vec![format_simple_card(&kw, &val_str, self.comment.as_deref())]
            }
            CardValue::Str(s) => format_string_cards(&kw, s, self.comment.as_deref()),
        }
    }
}

fn normalize_keyword(kw: &str) -> String {
    let trimmed = kw.trim();
    let mut out = String::with_capacity(trimmed.len().min(8));
    for c in trimmed.chars().take(8) {
        let uc = c.to_ascii_uppercase();
        if uc.is_ascii_alphanumeric() || uc == '-' || uc == '_' {
            out.push(uc);
        }
    }
    out
}

fn sanitize_ascii(s: &str) -> String {
    s.chars()
        .map(|c| {
            let b = c as u32;
            if (0x20..=0x7E).contains(&b) {
                c
            } else {
                ' '
            }
        })
        .collect()
}

fn format_fits_float(v: f64) -> String {
    if !v.is_finite() {
        return "0.0".to_string();
    }
    if v == 0.0 {
        return "0.0".to_string();
    }
    let abs = v.abs();
    if (1e-4..1e12).contains(&abs) {
        // Preserve up to 16 significant digits in fixed-point when it fits within 20 chars
        // (FITS fixed-format field width, preserving sub-microsecond precision on ~2e8 s timestamps).
        let int_digits = (abs.log10().floor() as i32 + 1).max(1);
        let frac = (16 - int_digits).clamp(1, 15) as usize;
        let mut s = format!("{:.prec$}", v, prec = frac);
        if s.contains('.') {
            while s.ends_with('0') {
                s.pop();
            }
            if s.ends_with('.') {
                s.push('0');
            }
        }
        if s.len() <= 20 {
            return s;
        }
    }
    // Scientific notation with 16 significant digits (matching CFITSIO %.15E),
    // uppercase E, and explicit exponent sign.
    let s = format!("{:.15E}", v);
    if let Some((mant, exp)) = s.split_once('E') {
        let mut m = mant.to_string();
        if m.contains('.') {
            while m.ends_with('0') {
                m.pop();
            }
            if m.ends_with('.') {
                m.push('0');
            }
        }
        let exp_num: i32 = exp.parse().unwrap_or(0);
        format!("{}E{:+03}", m, exp_num)
    } else {
        s
    }
}

fn format_simple_card(kw: &str, val_str: &str, comment: Option<&str>) -> [u8; 80] {
    let mut card = [b' '; 80];
    let kb = kw.as_bytes();
    card[..kb.len().min(8)].copy_from_slice(&kb[..kb.len().min(8)]);
    card[8] = b'=';
    card[9] = b' ';

    let vb = val_str.as_bytes();
    let vlen = vb.len().min(70);
    card[10..10 + vlen].copy_from_slice(&vb[..vlen]);

    if let Some(c) = comment {
        let clean = sanitize_ascii(c);
        if !clean.is_empty() {
            let start = (10 + vlen).max(30);
            if start + 3 < 80 {
                card[start] = b' ';
                card[start + 1] = b'/';
                card[start + 2] = b' ';
                let max_c = 80 - (start + 3);
                let cb = clean.as_bytes();
                let clen = cb.len().min(max_c);
                card[start + 3..start + 3 + clen].copy_from_slice(&cb[..clen]);
            }
        }
    }
    card
}

fn escape_fits_string(s: &str) -> String {
    let clean = sanitize_ascii(s);
    clean.replace('\'', "''")
}

fn format_string_cards(kw: &str, value: &str, comment: Option<&str>) -> Vec<[u8; 80]> {
    let escaped = escape_fits_string(value);
    // Max escaped content inside single quotes on a single 80-byte card:
    // bytes 0..10 = "KEYWORD = ", byte 10 = '\'', content, closing '\'' -> max 68 bytes.
    if escaped.len() <= 68 {
        let mut padded = escaped;
        // FITS fixed-format string convention: pad short strings to at least 8 chars inside quotes.
        if padded.len() < 8 {
            padded.push_str(&" ".repeat(8 - padded.len()));
        }
        let quoted = format!("'{}'", padded);
        let mut card = [b' '; 80];
        let kb = kw.as_bytes();
        card[..kb.len().min(8)].copy_from_slice(&kb[..kb.len().min(8)]);
        card[8] = b'=';
        card[9] = b' ';
        let qb = quoted.as_bytes();
        card[10..10 + qb.len()].copy_from_slice(qb);

        if let Some(c) = comment {
            let clean = sanitize_ascii(c);
            if !clean.is_empty() {
                let start = (10 + qb.len()).max(30);
                if start + 3 < 80 {
                    card[start] = b' ';
                    card[start + 1] = b'/';
                    card[start + 2] = b' ';
                    let max_c = 80 - (start + 3);
                    let cb = clean.as_bytes();
                    let clen = cb.len().min(max_c);
                    card[start + 3..start + 3 + clen].copy_from_slice(&cb[..clen]);
                }
            }
        }
        return vec![card];
    }

    // FITS 4.0 CONTINUE long-string convention.
    let mut cards = Vec::new();
    let bytes = escaped.as_bytes();
    let mut pos = 0;
    let mut first = true;
    while pos < bytes.len() {
        let remaining = bytes.len() - pos;
        let max_chunk = if remaining <= 68 { 68 } else { 67 };
        let mut end = (pos + max_chunk).min(bytes.len());
        // Do not split an escaped single-quote pair `''` across cards.
        if end < bytes.len() && bytes[end - 1] == b'\'' {
            let mut quote_run = 0usize;
            let mut k = end;
            while k > pos && bytes[k - 1] == b'\'' {
                quote_run += 1;
                k -= 1;
            }
            if quote_run % 2 == 1 && end - 1 > pos {
                end -= 1;
            }
        }
        let is_last = end >= bytes.len();
        let chunk = std::str::from_utf8(&bytes[pos..end]).unwrap_or("");
        let mut card = [b' '; 80];
        if first {
            let kb = kw.as_bytes();
            card[..kb.len().min(8)].copy_from_slice(&kb[..kb.len().min(8)]);
            card[8] = b'=';
            card[9] = b' ';
            first = false;
        } else {
            card[..8].copy_from_slice(b"CONTINUE");
            card[8] = b' ';
            card[9] = b' ';
        }
        let quoted = if is_last {
            format!("'{}'", chunk)
        } else {
            format!("'{}&'", chunk)
        };
        let qb = quoted.as_bytes();
        card[10..10 + qb.len()].copy_from_slice(qb);

        if is_last {
            if let Some(c) = comment {
                let clean = sanitize_ascii(c);
                if !clean.is_empty() {
                    let start = (10 + qb.len()).max(30);
                    if start + 3 < 80 {
                        card[start] = b' ';
                        card[start + 1] = b'/';
                        card[start + 2] = b' ';
                        let max_c = 80 - (start + 3);
                        let cb = clean.as_bytes();
                        let clen = cb.len().min(max_c);
                        card[start + 3..start + 3 + clen].copy_from_slice(&cb[..clen]);
                    }
                }
            }
        }
        cards.push(card);
        pos = end;
    }
    cards
}

/// A parsed or constructed FITS HDU header.
#[derive(Debug, Clone, Default)]
pub struct Header {
    cards: Vec<Card>,
}

impl Header {
    pub fn new() -> Self {
        Self { cards: Vec::new() }
    }

    pub fn with_cards(cards: Vec<Card>) -> Self {
        Self { cards }
    }

    /// Parse a FITS header starting at the beginning of `buf`.
    /// Returns `(Header, bytes_consumed)` where `bytes_consumed` is a multiple of 2880.
    pub fn parse(buf: &[u8]) -> Result<(Self, usize)> {
        if buf.len() < FITS_BLOCK_SIZE {
            return Err(FitsError::InvalidBlockAlignment(format!(
                "buffer length {} is smaller than one 2880-byte FITS block",
                buf.len()
            )));
        }

        let mut cards: Vec<Card> = Vec::new();
        let mut offset = 0usize;
        let mut found_end = false;

        while offset + 80 <= buf.len() {
            let raw = &buf[offset..offset + 80];
            offset += 80;

            // Check END card.
            if &raw[0..8] == b"END     " {
                found_end = true;
                break;
            }

            let kw = std::str::from_utf8(&raw[0..8])
                .map_err(|_| FitsError::InvalidHeader("non-UTF8 keyword".into()))?
                .trim_end()
                .to_string();

            // CONTINUE card appends to previous string card if that card ended with '&'.
            if kw == "CONTINUE" && &raw[8..10] == b"  " {
                if let Some(prev) = cards.last_mut() {
                    if let CardValue::Str(ref mut prev_s) = prev.value {
                        if prev_s.ends_with('&') {
                            let rest = std::str::from_utf8(&raw[10..80]).unwrap_or("");
                            if let Some((cont_str, cont_comment)) = parse_quoted_string(rest) {
                                prev_s.pop(); // remove trailing '&'
                                prev_s.push_str(&cont_str);
                                if let Some(c) = cont_comment {
                                    if prev.comment.is_none() {
                                        prev.comment = Some(c);
                                    }
                                }
                                continue;
                            }
                        }
                    }
                }
            }

            // Skip completely blank cards.
            if kw.is_empty() && raw[8..80].iter().all(|&b| b == b' ') {
                continue;
            }

            if kw == "COMMENT" || kw == "HISTORY" || &raw[8..10] != b"= " {
                let text = std::str::from_utf8(&raw[8..80])
                    .unwrap_or("")
                    .trim_end()
                    .to_string();
                cards.push(Card {
                    keyword: kw,
                    value: CardValue::Empty,
                    comment: if text.is_empty() { None } else { Some(text) },
                });
                continue;
            }

            let rest = std::str::from_utf8(&raw[10..80])
                .map_err(|_| FitsError::InvalidHeader(format!("non-UTF8 value for {kw}")))?;
            let (value, comment) = parse_value_and_comment(rest);
            cards.push(Card {
                keyword: kw,
                value,
                comment,
            });
        }

        if !found_end {
            return Err(FitsError::InvalidHeader(
                "missing END card in FITS header".into(),
            ));
        }

        // Trim trailing spaces on all parsed string values (FITS standard: trailing spaces in strings are insignificant).
        for card in &mut cards {
            if let CardValue::Str(ref mut s) = card.value {
                let trimmed_len = s.trim_end_matches(' ').len();
                s.truncate(trimmed_len);
            }
        }

        let blocks = offset.div_ceil(FITS_BLOCK_SIZE);
        let bytes_consumed = blocks * FITS_BLOCK_SIZE;
        if bytes_consumed > buf.len() {
            return Err(FitsError::InvalidBlockAlignment(format!(
                "header requires {bytes_consumed} bytes, buffer has {}",
                buf.len()
            )));
        }

        Ok((Self { cards }, bytes_consumed))
    }

    /// Serialize this header into a 2880-byte-aligned FITS header byte vector ending with `END`.
    pub fn to_fits_bytes(&self) -> Vec<u8> {
        let has_longstrn = self.get("LONGSTRN").is_some();
        let mut uses_continue = false;
        let mut phys_cards: Vec<[u8; 80]> = Vec::with_capacity(self.cards.len() + 2);
        for card in &self.cards {
            if card.keyword.eq_ignore_ascii_case("END") {
                continue;
            }
            let rendered = card.to_fits_cards();
            if rendered.len() > 1
                && !card.keyword.eq_ignore_ascii_case("COMMENT")
                && !card.keyword.eq_ignore_ascii_case("HISTORY")
            {
                uses_continue = true;
            }
            phys_cards.extend(rendered);
        }
        if uses_continue && !has_longstrn {
            let ls = Card::new_str(
                "LONGSTRN",
                "OGIP 1.0",
                Some("The OGIP long string convention may be used".to_string()),
            );
            // Insert before CHECKSUM/DATASUM if present, else before END.
            let insert_pos = phys_cards
                .iter()
                .position(|c| &c[0..8] == b"DATASUM " || &c[0..8] == b"CHECKSUM")
                .unwrap_or(phys_cards.len());
            for (k, c) in ls.to_fits_cards().into_iter().enumerate() {
                phys_cards.insert(insert_pos + k, c);
            }
        }
        let mut end_card = [b' '; 80];
        end_card[0..3].copy_from_slice(b"END");
        phys_cards.push(end_card);

        let raw_len = phys_cards.len() * 80;
        let padded_len = raw_len.div_ceil(FITS_BLOCK_SIZE) * FITS_BLOCK_SIZE;
        let mut out = vec![b' '; padded_len];
        for (i, c) in phys_cards.iter().enumerate() {
            out[i * 80..(i + 1) * 80].copy_from_slice(c);
        }
        out
    }

    pub fn cards(&self) -> &[Card] {
        &self.cards
    }

    pub fn cards_mut(&mut self) -> &mut Vec<Card> {
        &mut self.cards
    }

    pub fn get(&self, key: &str) -> Option<&Card> {
        self.cards
            .iter()
            .find(|c| c.keyword.eq_ignore_ascii_case(key))
    }

    pub fn get_str(&self, key: &str) -> Option<&str> {
        match &self.get(key)?.value {
            CardValue::Str(s) => Some(s.as_str()),
            _ => None,
        }
    }

    pub fn get_i64(&self, key: &str) -> Option<i64> {
        match &self.get(key)?.value {
            CardValue::Integer(n) => Some(*n),
            CardValue::Float(f)
                if f.is_finite()
                    && f.fract() == 0.0
                    && *f >= (i64::MIN as f64)
                    && *f < 9_223_372_036_854_775_808.0 =>
            {
                Some(*f as i64)
            }
            _ => None,
        }
    }

    pub fn get_f64(&self, key: &str) -> Option<f64> {
        match &self.get(key)?.value {
            CardValue::Float(f) => Some(*f),
            CardValue::Integer(n) => Some(*n as f64),
            _ => None,
        }
    }

    pub fn get_bool(&self, key: &str) -> Option<bool> {
        match &self.get(key)?.value {
            CardValue::Logical(b) => Some(*b),
            _ => None,
        }
    }

    pub fn extname(&self) -> Option<&str> {
        self.get_str("EXTNAME")
    }

    pub fn extver(&self) -> Option<i64> {
        self.get_i64("EXTVER")
    }

    /// Insert or update a card. `COMMENT` and `HISTORY` are always appended.
    pub fn set(&mut self, card: Card) {
        if card.keyword.eq_ignore_ascii_case("COMMENT")
            || card.keyword.eq_ignore_ascii_case("HISTORY")
        {
            self.cards.push(card);
            return;
        }
        if let Some(existing) = self
            .cards
            .iter_mut()
            .find(|c| c.keyword.eq_ignore_ascii_case(&card.keyword))
        {
            *existing = card;
        } else {
            self.cards.push(card);
        }
    }

    pub fn set_str(&mut self, key: &str, val: impl AsRef<str>, comment: &str) {
        let c = if comment.is_empty() {
            None
        } else {
            Some(comment.to_string())
        };
        self.set(Card::new_str(key, val.as_ref().to_string(), c));
    }

    pub fn set_i64(&mut self, key: &str, val: i64, comment: &str) {
        let c = if comment.is_empty() {
            None
        } else {
            Some(comment.to_string())
        };
        self.set(Card::new_int(key, val, c));
    }

    pub fn set_f64(&mut self, key: &str, val: f64, comment: &str) {
        let c = if comment.is_empty() {
            None
        } else {
            Some(comment.to_string())
        };
        self.set(Card::new_float(key, val, c));
    }

    pub fn set_bool(&mut self, key: &str, val: bool, comment: &str) {
        let c = if comment.is_empty() {
            None
        } else {
            Some(comment.to_string())
        };
        self.set(Card::new_bool(key, val, c));
    }

    pub fn add_comment(&mut self, text: impl Into<String>) {
        self.cards.push(Card::comment(text));
    }

    pub fn add_history(&mut self, text: impl Into<String>) {
        self.cards.push(Card::history(text));
    }

    pub fn remove(&mut self, key: &str) {
        self.cards.retain(|c| !c.keyword.eq_ignore_ascii_case(key));
    }
}

fn parse_quoted_string(rest: &str) -> Option<(String, Option<String>)> {
    let trimmed = rest.trim_start();
    let inner = trimmed.strip_prefix('\'')?;
    let mut out = String::new();
    let mut chars = inner.char_indices().peekable();
    while let Some((idx, ch)) = chars.next() {
        if ch == '\'' {
            if let Some(&(_, '\'')) = chars.peek() {
                chars.next();
                out.push('\'');
            } else {
                // Closing quote found.
                let after = &inner[idx + 1..];
                let comment = after.find('/').map(|pos| after[pos + 1..].trim().to_string());
                return Some((out, comment.filter(|c| !c.is_empty())));
            }
        } else {
            out.push(ch);
        }
    }
    // Unterminated quote: return what we have gracefully.
    Some((out, None))
}

fn parse_value_and_comment(rest: &str) -> (CardValue, Option<String>) {
    let trimmed = rest.trim_start();
    if trimmed.is_empty() {
        return (CardValue::Empty, None);
    }

    if trimmed.starts_with('\'') {
        if let Some((s, c)) = parse_quoted_string(trimmed) {
            return (CardValue::Str(s), c);
        }
    }

    let (val_part, comment) = match trimmed.find('/') {
        Some(pos) => {
            let c = trimmed[pos + 1..].trim();
            (
                trimmed[..pos].trim(),
                if c.is_empty() {
                    None
                } else {
                    Some(c.to_string())
                },
            )
        }
        None => (trimmed.trim(), None),
    };

    if val_part.is_empty() {
        return (CardValue::Empty, comment);
    }
    if val_part == "T" {
        return (CardValue::Logical(true), comment);
    }
    if val_part == "F" {
        return (CardValue::Logical(false), comment);
    }
    if let Ok(n) = val_part.parse::<i64>() {
        return (CardValue::Integer(n), comment);
    }
    // Handle Fortran 'D' exponent notation occasionally found in legacy FITS files.
    let norm_float = val_part.replace(['D', 'd'], "E");
    if let Ok(f) = norm_float.parse::<f64>() {
        return (CardValue::Float(f), comment);
    }

    (CardValue::Str(val_part.to_string()), comment)
}
