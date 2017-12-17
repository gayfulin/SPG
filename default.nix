with import <nixpkgs> {};
stdenv.mkDerivation {
  name = "SPG";
  src = ./.;
  buildInputs = [
    cmake
    gfortran
    valgrind
  ];
  separateDebugInfo = true;
  hardeningDisable = [ "all" ];
  doCheck = true;
}
