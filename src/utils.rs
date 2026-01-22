
use std::io::{Write};
use env_logger::{Builder};
use log::LevelFilter;
use std::{thread};
use sysinfo::{System, Pid, ProcessesToUpdate};

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

pub fn mem_logging() {
    let pid: u32 = std::process::id();

    // Spawn a thread to monitor memory
    thread::spawn(move || {
        use std::fs::File;
        use std::io::{BufWriter, Write};
        use std::time::{Duration, Instant};
        use log::info;

        let start = Instant::now();

        let file = File::create("memory.log")
            .expect("failed to create memory log");

        let mut stdout_mem = BufWriter::new(file);
        let mut system = System::new_all();

        loop {
            system.refresh_processes(
                ProcessesToUpdate::Some(&[Pid::from(pid as usize)]), 
                true
            );
            let elapsed = start.elapsed().as_secs();
            if let Some(process) = system.process(Pid::from(pid as usize)) {
                writeln!(stdout_mem,"{:.2}, {:.3} GB", elapsed, process.memory() as f32 / ( 1024 * 1024 * 1024 ) as f32 ).ok();
                info!("{:.2}, {:.3} GB", elapsed, process.memory() as f32 / ( 1024 * 1024 * 1024 ) as f32 );

            } else {
                writeln!(stdout_mem,"Process not found!").ok();
            }
            stdout_mem.flush().ok();
            thread::sleep(Duration::from_secs(1));
        }
    });
}