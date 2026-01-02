/// Returns the maximum resident set size (RSS) memory usage of the current process in megabytes.
///
/// This function uses the `getrusage` system call to retrieve resource usage information.
/// The unit of `ru_maxrss` varies by OS (kilobytes on Linux, bytes on macOS), so it adjusts
/// the calculation accordingly.
///
/// # Returns
/// The maximum memory usage in megabytes (`f64`).
///
/// # Safety
/// This function uses `unsafe` block to call the `libc::getrusage` C function.
/// It assumes the `getrusage` call is safe and the `rusage` struct is initialized correctly.
///
/// # Example
/// ```rust, ignore
/// // This example will show memory usage, but its value depends on the runtime environment.
/// use gxf2bed::max_mem_usage_mb; //
///
/// let mem_usage = max_mem_usage_mb();
/// println!("Max memory usage: {:.2} MB", mem_usage);
/// ```
pub fn max_mem_usage_mb() -> f64 {
    let rusage = unsafe {
        let mut rusage = std::mem::MaybeUninit::uninit();
        libc::getrusage(libc::RUSAGE_SELF, rusage.as_mut_ptr());
        rusage.assume_init()
    };
    let maxrss = rusage.ru_maxrss as f64;
    if cfg!(target_os = "macos") {
        maxrss / 1024.0 / 1024.0
    } else {
        maxrss / 1024.0
    }
}
