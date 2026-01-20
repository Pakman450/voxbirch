
use std::io::{Write};
use env_logger::{Builder};
use log::LevelFilter;

pub fn calc_time_breakdown (duration_mark: &std::time::Duration) -> (
    u64,
    u64,
    u64,
    u64,
    u32
) {

    let total_secs = duration_mark.as_secs();
    let hours = total_secs / 3600;
    let minutes = (total_secs % 3600) / 60;
    let seconds = total_secs % 60;
    let milliseconds = duration_mark.subsec_millis();

    // If the process takes less than 1 sec, 
    // you have to pass milliseconds as a
    // total duration, which is used for calculating
    // mol proccessed per sec
    if total_secs == 0 {

        return  (
            milliseconds as u64,
            hours,
            minutes,
            seconds,
            milliseconds
        )   
    }


    (
        total_secs,
        hours,
        minutes,
        seconds,
        milliseconds
    )
}

pub fn init_logging(verbosity: u8) {
    let mut builder = Builder::new();

    match verbosity {
        0 => builder.filter_level(LevelFilter::Warn),   // default
        1 => builder.filter_level(LevelFilter::Info),
        2 => builder.filter_level(LevelFilter::Debug),
        _ => builder.filter_level(LevelFilter::Trace),
    };

    builder.format(|buf, record| {

        let level_style = buf.default_level_style(record.level());
        let level = level_style.value(record.level());

        let file = record.file().unwrap_or("unknown");
        let line = record.line().unwrap_or(0);

        writeln!(
            buf,
            "[{} {}:{} {}] {}",
            level,
            file,
            line,
            record.target(),
            record.args()
        )
    });
    
    builder.init();
}