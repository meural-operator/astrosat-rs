//! TOML/JSON configuration parsing and normalization endpoint (`POST /api/config/parse`).
//!
//! Uses `astrosat_core::toml_lite` so UI config imports share 100% identical TOML semantics
//! and `AppConfig` schema validation with `astrosat --config`.

use std::io::Read;

use astrosat_core::{parse_toml_to_json, to_canonical_toml, AppConfig};
use serde::{Deserialize, Serialize};
use serde_json::Value;
use tiny_http::Request;

use crate::http_util::{ApiError, ApiResult, MAX_BODY_BYTES};

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct ConfigParseRequest {
    toml: Option<String>,
    content: Option<String>,
    format: Option<String>,
}

#[derive(Debug, Clone, Serialize)]
pub struct ConfigParseResponse {
    pub config: AppConfig,
    pub sections: Vec<String>,
    pub canonical_toml: String,
}

fn read_raw_body(req: &mut Request) -> ApiResult<String> {
    let mut buf = Vec::new();
    let reader = req.as_reader();
    let mut limited = reader.take((MAX_BODY_BYTES + 1) as u64);
    limited
        .read_to_end(&mut buf)
        .map_err(|e| ApiError::bad_request(format!("failed to read request body: {e}")))?;
    if buf.len() > MAX_BODY_BYTES {
        return Err(ApiError::bad_request(format!(
            "request body exceeds {} bytes limit",
            MAX_BODY_BYTES
        )));
    }
    if buf.is_empty() {
        return Err(ApiError::bad_request("empty request body"));
    }
    String::from_utf8(buf).map_err(|_| ApiError::bad_request("request body must be valid UTF-8"))
}

fn content_type_media(req: &Request) -> Option<String> {
    req.headers()
        .iter()
        .find(|h| h.field.equiv("Content-Type"))
        .map(|h| {
            h.value
                .as_str()
                .split(';')
                .next()
                .unwrap_or("")
                .trim()
                .to_ascii_lowercase()
        })
}

/// Handle `POST /api/config/parse`.
pub fn handle_config_parse(req: &mut Request) -> ApiResult<ConfigParseResponse> {
    let raw = read_raw_body(req)?;
    let ct = content_type_media(req);

    let val: Value = match ct.as_deref() {
        Some("text/plain") | Some("application/toml") => parse_toml_to_json(&raw)
            .map_err(|e| ApiError::bad_request(format!("invalid TOML config: {e}")))?,
        _ => {
            let parsed: ConfigParseRequest = serde_json::from_str(&raw)
                .map_err(|e| ApiError::bad_request(format!("invalid JSON body: {e}")))?;
            match (parsed.toml, parsed.content) {
                (Some(_), Some(_)) => {
                    return Err(ApiError::bad_request(
                        "provide either 'toml' or 'content', not both",
                    ));
                }
                (Some(s), None) => {
                    if let Some(ref f) = parsed.format {
                        if !f.trim().eq_ignore_ascii_case("toml") {
                            return Err(ApiError::bad_request(format!(
                                "field 'toml' cannot be combined with format '{}'",
                                f.trim()
                            )));
                        }
                    }
                    parse_toml_to_json(&s)
                        .map_err(|e| ApiError::bad_request(format!("invalid TOML config: {e}")))?
                }
                (None, Some(s)) => {
                    let fmt = parsed
                        .format
                        .as_deref()
                        .unwrap_or("toml")
                        .trim()
                        .to_ascii_lowercase();
                    match fmt.as_str() {
                        "toml" => parse_toml_to_json(&s).map_err(|e| {
                            ApiError::bad_request(format!("invalid TOML config: {e}"))
                        })?,
                        "json" => serde_json::from_str::<Value>(
                            s.strip_prefix('\u{FEFF}').unwrap_or(&s),
                        )
                        .map_err(|e| {
                            ApiError::bad_request(format!("invalid JSON config: {e}"))
                        })?,
                        "auto" => {
                            let s_clean = s.strip_prefix('\u{FEFF}').unwrap_or(&s);
                            if s_clean.trim_start().starts_with('{') {
                                serde_json::from_str::<Value>(s_clean).map_err(|e| {
                                    ApiError::bad_request(format!("invalid JSON config: {e}"))
                                })?
                            } else {
                                parse_toml_to_json(s_clean).map_err(|e| {
                                    ApiError::bad_request(format!("invalid TOML config: {e}"))
                                })?
                            }
                        }
                        other => {
                            return Err(ApiError::bad_request(format!(
                                "unsupported config format '{other}' (expected 'toml' or 'json')"
                            )));
                        }
                    }
                }
                (None, None) => {
                    return Err(ApiError::bad_request(
                        "missing 'toml' or 'content' field in request body",
                    ));
                }
            }
        }
    };

    let Some(map) = val.as_object() else {
        return Err(ApiError::bad_request("config root must be a table/object"));
    };

    let mut sections = Vec::new();
    for name in ["czti", "laxpc", "bary", "timing", "ui"] {
        if map.contains_key(name) {
            sections.push(name.to_string());
        }
    }

    let cfg: AppConfig = serde_json::from_value(val)
        .map_err(|e| ApiError::bad_request(format!("config schema error: {e}")))?;

    cfg.validate()
        .map_err(|e| ApiError::bad_request(format!("config validation error: {e}")))?;

    let canonical_toml = to_canonical_toml(&cfg);

    Ok(ConfigParseResponse {
        config: cfg,
        sections,
        canonical_toml,
    })
}
