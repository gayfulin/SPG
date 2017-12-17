with import <nixpkgs> {};
stdenv.mkDerivation {
  name = "SPG";
  src = ./.;
  buildInputs = [
    gfortran
  ];
  buildPhase = ''
    ./gcompile
  '';
  hardeningDisable = "all";
  installPhase = ''
    mkdir -p $out/bin
    cp main.exe $out/bin/SPG.exe
  '';
}
