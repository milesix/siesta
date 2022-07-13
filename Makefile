default: siesta

siesta:
	(cd Src; make)

utils:
	(cd Util/Denchar/Src; make)
	(cd Util/VCA; make)
	(cd Util/Grid; make)
	(cd Util/Gen-basis; make)
	(cd Util/COOP; make)
	(cd Util/TS/TBtrans; make)
	(cd Util/TS/tscontour; make)
	(cd Util/TS/tshs2tshs; make)
	(cd Util/TS/ts2ts; make)
	(cd Util/Vibra/Src; make)
	(cd Util/Macroave/Src; make)
	(cd Util/STM/ol-stm/Src; make)
	(cd Util/Unfolding/Src; make)
	(cd Pseudo/converters/psml2psf; make)

clean_utils:
	(cd Util/Denchar/Src; make clean)
	(cd Util/VCA; make clean)
	(cd Util/Grid; make clean)
	(cd Util/Gen-basis; make clean)
	(cd Util/COOP; make clean)
	(cd Util/TS/TBtrans; make clean)
	(cd Util/TS/tscontour; make clean)
	(cd Util/TS/tshs2tshs; make clean)
	(cd Util/TS/ts2ts; make clean)
	(cd Util/Vibra/Src; make clean)
	(cd Util/Macroave/Src; make clean)
	(cd Util/STM/ol-stm/Src; make clean)
	(cd Util/Unfolding/Src; make clean)
	(cd Pseudo/converters/psml2psf; make clean)

clean_siesta:
	(cd Src; make clean)

clean: clean_utils clean_siesta

