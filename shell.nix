{ pkgs ? import <nixpkgs> {} } : with pkgs; mkShell {
	packages = [
		(python3.withPackages (ps: with ps; [
			numpy
			matplotlib
			qutip
			scipy
			tqdm
		]))
	];
}
