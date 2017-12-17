with import <nixpkgs> {};
stdenv.mkDerivation {
  name = "SPG";
  src = ./.;
  buildInputs = [
    cmake
    gfortran
  ];
  hardeningDisable = [ "all" ];
}
