# RareXSec

## Build instructions

1. Clone the repository and enter it:
   ```bash
   git clone https://github.com/<your-org>/rarexsec.git
   cd rarexsec
   ```
2. Build the project libraries from the provided `build/Makefile`:
   ```bash
   cd build
   make -j"$(nproc)"
   ```
   The compiled artifacts are placed under `build/lib` and `build/bin`.  Add
   `-j"$(nproc)"` if you want to compile in parallel.
3. (Optional) Install the headers, libraries, and helper scripts into a prefix:
   ```bash
   make install PREFIX=/desired/install/location
   ```

## Running the example ROOT macro

After building, the rarexsec libraries must be discoverable by ROOT.  The `rarexsec-root.sh` wrapper sets up the include and library paths and executes any macro you pass to it.  From the repository root run:

```bash
./rarexsec-root.sh -b -q macros/example_macro.C
```

This command loads the generated libraries, configures ROOT include paths, and runs the `example_macro()` entry point defined in `macros/example_macro.C`.

Alternatively, you can call the setup macro manually from ROOT by providing the absolute library and include directories:

```bash
root -l -q -e 'setup_rarexsec("$PWD/build/lib/librarexsec.so","$PWD/include")' macros/example_macro.C
```