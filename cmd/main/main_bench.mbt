// main_bench.mbt

// Simple test that measures runtime using `now_seconds()` and prints it.
test {
  let start = now_seconds()
  run_simple_simulation() |> ignore
  let end = now_seconds()
  if start > 0.0 && end > start {
    let elapsed = end - start
    println("Simulation runtime (seconds): \{elapsed}")
  } else {
    println("Simulation runtime: timing API unavailable")
  }
}
