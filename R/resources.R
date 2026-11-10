#' Threads and memory a plink4r call may use
#'
#' Read-only detection, used whenever a function is called with `threads = 0` (or `mem_gb = 0` where it takes
#' one). It covers the common setups only. Threads: the scheduler's allocation (`SLURM_CPUS_PER_TASK`,
#' `SLURM_CPUS_ON_NODE`; PBS/Torque `NCPUS`, `PBS_NP`; LSF `LSB_DJOB_NUMPROC`; SGE `NSLOTS`), capped by the
#' CPU affinity mask; with no scheduler, the affinity mask when it is narrower than the machine (a cgroup or
#' taskset limit), else 4. The total core count of a node is never used. Windows and macOS: 4. Memory: half of
#' the smallest of the SLURM allocation (`SLURM_MEM_PER_NODE`, or `SLURM_MEM_PER_CPU` times the threads), the
#' cgroup limit and `MemAvailable`, at most 16 GB on a machine with neither (a shared workstation or login
#' node); 2 GB when nothing is known. Nothing is written: no environment variable, no global thread setting.
#'
#' @param threads Threads; 0 detects, a positive value is returned as given.
#' @param mem_gb Memory budget in GB; 0 detects, a positive value is returned as given.
#' @return A list with `threads`, `mem_gb`, and where each came from (`threads_from`, `mem_from`).
#' @export
plink4r_resources <- function(threads = 0L, mem_gb = 0) {
  if (!is.numeric(threads) || length(threads) != 1L || !is.finite(threads) || threads < 0 || threads != round(threads))
    stop("threads must be 0 (detect) or a positive integer.", call. = FALSE)
  if (!is.numeric(mem_gb) || length(mem_gb) != 1L || !is.finite(mem_gb) || mem_gb < 0)
    stop("mem_gb must be 0 (detect) or a positive number.", call. = FALSE)
  resources_cpp(as.integer(threads), as.numeric(mem_gb))
}

# threads for a call: 0 resolves through plink4r_resources(), a positive integer is kept.
.threads <- function(threads) as.integer(plink4r_resources(threads, 1)$threads)
# memory budget (GB) for a call: 0 resolves through plink4r_resources() given the resolved threads.
.mem_gb <- function(mem_gb, threads) as.numeric(plink4r_resources(threads, mem_gb)$mem_gb)
