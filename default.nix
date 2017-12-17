with import <nixpkgs> {};
stdenv.mkDerivation {
  name = "SPG";
  src = ./.;
  buildInputs = [
    cmake
    gfortran
    valgrind
  ];
  separateDebugInfo = stdenv.isLinux; # Currently unsupported on darwin
  hardeningDisable = [ "all" ];
  doCheck = true;
}
