
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