# Performance benchmarks

Some notable function in Pymech are benchmarked across different versions.
The benchmarks are run using the tool [airspeed velocity /
`asv`](https://asv.readthedocs.io).


```{note}
The results of the benchmark are <a href="index.html">available here</a>
```

##  How to run benchmarks

- Either install `pymech` with `[dev]` extras or `pip install asv pyperf virtualenv`.

- Change to the directory `docs/asv_bench`

- Close all heavy applications (browsers, video players, chat clients etc.) and
  tune the system for reducing jitter during benchmarking:
  ```
  sudo python -m pyperf system tune
  ```

- Execute the benchmarks
	- on the current commit with the default Python:

          asv run

	- on the current commit with the PyPy installed using conda

          asv run --environment 'conda:*=*pypy'

	- on all git tags:
      ```bash
      ./asv_run_all_tags.sh
      ./asv_run_all_tags.sh --environment 'conda:*=*pypy'
      ```

- Build HTML of all the results `asv publish`

- View the results `asv preview`
