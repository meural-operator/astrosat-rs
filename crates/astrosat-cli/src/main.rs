//! Unified `astrosat` CLI entrypoint.

mod bary_cmd;
mod cli;
mod config;
mod czti_cmd;
mod inspect;
mod laxpc_cmd;
mod output;
mod timing;
mod validate;

use std::process::ExitCode;

use anyhow::Context;
use clap::{ArgMatches, CommandFactory, FromArgMatches};
use tracing_subscriber::EnvFilter;

use crate::bary_cmd::{print_bary_human, run_bary};
use crate::cli::{is_cli_arg, leaf_matches, resolve_flag, Cli, TopCommand};
use crate::config::{load_app_config, AppConfig};
use crate::czti_cmd::dispatch_czti;
use crate::inspect::{print_inspect_human, run_inspect};
use crate::laxpc_cmd::dispatch_laxpc;
use crate::output::{format_error_hint, OutputMode};
use crate::timing::{print_timing_human, run_timing, TimingParams};
use crate::validate::{print_validate_human, run_validate};

fn main() -> ExitCode {
    let matches = Cli::command().get_matches();
    let cli = match Cli::from_arg_matches(&matches) {
        Ok(c) => c,
        Err(e) => e.exit(),
    };
    let json_mode = cli.json;
    init_tracing(cli.verbose, cli.quiet, json_mode);

    match run(cli, &matches) {
        Ok(exit_code) => exit_code,
        Err(err) => {
            if json_mode {
                let payload = serde_json::json!({
                    "error": format!("{err:#}"),
                    "hint": format_error_hint(&err),
                });
                eprintln!("{}", serde_json::to_string_pretty(&payload).unwrap_or_default());
            } else {
                eprintln!("error: {err:#}");
                if let Some(hint) = format_error_hint(&err) {
                    eprintln!("{hint}");
                }
            }
            ExitCode::FAILURE
        }
    }
}

fn init_tracing(verbose: bool, quiet: bool, json: bool) {
    let level = if quiet || json {
        "warn"
    } else if verbose {
        "debug"
    } else {
        "info"
    };
    let filter = EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new(level));
    let _ = tracing_subscriber::fmt()
        .with_env_filter(filter)
        .with_writer(std::io::stderr)
        .with_target(false)
        .try_init();
}

fn run(cli: Cli, root_matches: &ArgMatches) -> anyhow::Result<ExitCode> {
    let leaf = leaf_matches(root_matches);
    let app_cfg = if let Some(ref cfg_path) = cli.config {
        load_app_config(cfg_path)?
    } else {
        AppConfig::default()
    };

    if let Some(threads) = cli.threads.or(app_cfg.threads) {
        if threads == 0 {
            anyhow::bail!("--threads must be >= 1");
        }
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .context("failed to configure Rayon global thread pool")?;
    }

    let out = OutputMode {
        json: cli.json,
        quiet: cli.quiet,
    };

    match &cli.command {
        TopCommand::Czti { command } => {
            dispatch_czti(command, &app_cfg.czti, leaf, out)?;
        }
        TopCommand::Laxpc { command } => {
            dispatch_laxpc(command, &app_cfg.laxpc, leaf, out)?;
        }
        TopCommand::Bary(args) => {
            let report = run_bary(args, &app_cfg.bary, leaf)?;
            out.emit(&report, print_bary_human)?;
        }
        TopCommand::Inspect(args) => {
            let report = run_inspect(&args.files, args.full)?;
            out.emit(&report, print_inspect_human)?;
        }
        TopCommand::Validate(args) => {
            let report = run_validate(&args.files, args.strict);
            out.emit(&report, print_validate_human)?;
            if !report.passed {
                return Ok(ExitCode::from(2));
            }
        }
        TopCommand::Timing(args) => {
            let params = TimingParams {
                input: args.input.clone(),
                hdu: args.hdu.clone(),
                gti_file: args.gti.clone(),
                freq_hz: if args.freq.is_some() {
                    args.freq
                } else if args.period.is_some() {
                    None
                } else {
                    app_cfg.timing.freq_hz
                },
                period_sec: args.period,
                fdot: if is_cli_arg(leaf, "fdot") {
                    args.fdot
                } else {
                    app_cfg.timing.fdot.unwrap_or(args.fdot)
                },
                fddot: if is_cli_arg(leaf, "fddot") {
                    args.fddot
                } else {
                    app_cfg.timing.fddot.unwrap_or(args.fddot)
                },
                epoch: args.epoch,
                freq_min: args.freq_min,
                freq_max: args.freq_max,
                search_span_hz: args.search_span,
                df_hz: args.df,
                n_trials: args.n_trials,
                harmonics: if is_cli_arg(leaf, "harmonics") {
                    args.harmonics
                } else {
                    app_cfg.timing.harmonics.unwrap_or(args.harmonics)
                },
                nbins: if is_cli_arg(leaf, "nbins") {
                    args.nbins
                } else {
                    app_cfg.timing.nbins.unwrap_or(args.nbins)
                },
                pi_min: args.pi_min,
                pi_max: args.pi_max,
                unit_filter: args.unit,
                layer_filter: args.layer,
            };
            let res = run_timing(&params)?;
            out.emit(&res, print_timing_human)?;
        }
        TopCommand::Ui(args) => {
            let host = if is_cli_arg(leaf, "host") {
                args.host.clone()
            } else {
                app_cfg.ui.host.clone().unwrap_or_else(|| args.host.clone())
            };
            let port = if is_cli_arg(leaf, "port") {
                args.port
            } else {
                app_cfg.ui.port.unwrap_or(args.port)
            };
            let data_dir = if let Some(dir) = args.data_dir.clone().or_else(|| app_cfg.ui.data_dir.clone()) {
                if !dir.is_dir() {
                    anyhow::bail!("--data-dir is not an existing directory: {}", dir.display());
                }
                dir
            } else {
                std::env::current_dir().unwrap_or_else(|_| std::path::PathBuf::from("."))
            };
            let open_browser = resolve_flag(leaf, "open", "no_open", app_cfg.ui.open_browser, true);
            astrosat_ui::serve(astrosat_ui::UiConfig {
                host,
                port,
                data_dir,
                open_browser,
                threads: cli.threads.or(app_cfg.threads),
            })?;
        }
    }

    Ok(ExitCode::SUCCESS)
}
