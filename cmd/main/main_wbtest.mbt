// benchmark_test.mbt
// Lightweight smoke test for the simulation entrypoint.
fn wrapped_main() -> Int {
  run_simple_simulation() |> ignore
  0
}

test {
  // Execute the wrapped simple simulation to ensure it runs without panics.
  wrapped_main() |> ignore
}