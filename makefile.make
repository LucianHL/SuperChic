LIBFLAGlha = -lLHAPDF
LHAPDFLIB = `lhapdf-config --prefix`/lib
FC = gfortran -g -w -Wall -fbacktrace -fcheck=all -fPIC
# -ffixed-line-length-132

#####################

HOME = $(PWD)
SOURCEDIR = $(PWD)/src
INCPATH = $(SOURCEDIR)/inc

OBJ_PATH = $(PWD)/obj/

$(shell mkdir -p lib/)

# FFLAGS 	= -fno-automatic -fno-f2c -O2 -g  -I$(INCPATH) 
FFLAGS 	= -fno-automatic -fno-f2c -fPIC -O2 -g  -I$(INCPATH) -fimplicit-none -Wall -Wextra -Wno-unused-dummy-argument

DIRS	 =	$(SOURCEDIR)/int:\
		$(SOURCEDIR)/main:\
		$(SOURCEDIR)/mes:\
		$(SOURCEDIR)/PDFs:\
		$(SOURCEDIR)/phase:\
		$(SOURCEDIR)/sPDFs:\
		$(SOURCEDIR)/subamps:\
		$(SOURCEDIR)/surv:\
		$(SOURCEDIR)/EW:\
		$(SOURCEDIR)/unw:\
		$(SOURCEDIR)/user:\
		$(SOURCEDIR)/int:\
		$(SOURCEDIR)/var:\
		$(SOURCEDIR)/init:\
		$(SOURCEDIR)/LbyL:\
		$(SOURCEDIR)/inition:\
		$(SOURCEDIR)/diss:\
		$(SOURCEDIR)/DHELAS:\
		$(SOURCEDIR)/MODEL:\
		$(SOURCEDIR)/MG:\

VPATH = $(DIRS)

#############

Dissf = \
CLAS.o \
gd-fit-11.o \
CB.o \
Elastic.o \
F1F2.o \
F2ev.o \
rfit.o \
alphaEM.o \
qinit.o \
apfelinit.o \

LbyLf = \
B0F2M.o \
D0404M.o \
xspenz.o \
eett_aux.o \
C01_gen.o \

Mesf = \
calcmes.o \
mesint.o \
wfinit.o \
wfoctet.o \
wfsinglet.o \

PDFsfLHA = \
alphas.o \
inpdf.o \

PDFsfUSER = \
alphasuser.o \
inpdfuser.o \

Intf = \
vegas.o \
rann.o \

Phasef = \
2body.o \
2bodyw.o \
2jetps.o \
2jetpsm.o \
3jetps.o \
boost.o \
chic0decay3.o \
chic1decay3.o \
chic1decay2s.o \
chic1decay2f.o \
chic2decay3.o \
chic2decay2s.o \
chic2decay2f.o \
jpsidecayphot.o \
genpol1.o \
genpol1rf.o \
genpol2.o \
rambo.o \
6body.o \
6bodyinit.o \
4body.o \
4bodyinit.o \
3body.o \
3bodyinit.o \
2bodyinit.o \
wwcorr.o \
jpsidecay.o \
rhodecay.o \
chidecay.o \
monow.o \
alpdecay.o \
pAboost.o \
pAinit.o \
AAinit.o \
genpol_axial.o \

Subampsf = \
chi0.o \
chi1.o \
chi2.o \
etaq.o \
higgs.o \
higgsinit.o \
pipi.o \
qqjets.o \
diphoton.o \
etaeta.o \
etapetap.o \
etaetap.o \
eta.o \
gggjets.o \
pipixy.o \
rhorho.o \
djpsi.o \
djpsip.o \
djpsipp.o \
ggjets.o \
qqgjets.o \
rhorhoxy.o \
wwpol.o \
llpol.o \
mhv.o \
lightlightpol.o \
higgsgam.o \
higgsgaminit.o \
alp.o \
mmpol.o \
monop.o \
lloff.o \
spinors.o \
wwoff.o \
wwoff_amp.o \

Survf = \
initparsr.o \
formfac.o \
formfacphot.o \
formfacgam.o \
seik.o \
seikphot.o\
seikgam.o \
screeningint.o \
readscreen.o \
formfacgamel.o \
formfacgamion.o \
formfacgamionp.o \
tpint.o \
seikgamion.o \
screeningionint.o \
seikion.o \
betaionex.o \
betaion.o \
seikphotionp.o \
formfacphotionp.o \
formfacgamoff.o \
formfacgamoff_ion.o \
formfacgamoff_ion_surv.o \
formfacgamoff_surv.o \
s2_diss_int.o \
s2_diss_init.o \

Userf = \
cuts.o \
histo.o \

Mainf = \
bare.o \
header.o \
main.o \
process.o \
wtgen.o \
wtgengam.o \
header_out.o \

sPDFsf = \
calchg.o \
hpdfint.o \
calcsud.o \
sudint.o \
sPDF.o \

Unwf = \
unweight.o \
headerlhe.o \
unweight_q.o \
unwprint_q.o \
unwprint.o \
unweight.o \

Varf = \
mu.o \
nf.o \
string.o \
varfuncs.o \
gaminit.o \
gauss.o \

InitfLHA = \
alphas.o \
initsud.o \
nf.o \
string.o \
hg.o \
inithg.o \
initpars.o \
calcop.o \
calcscreen.o \
opacityint.o \
screeningint.o \
screening.o \
opacity.o \
PDF.o \
PDFlha.o \
Sudakov.o \
inpdf.o \
dd.o \
Elastic.o \
sdcoh.o \
rann.o \
readscreen.o \
sdincoh.o \
F1F2.o \
F2ev.o \
rfit.o \
gd-fit-11.o \
CLAS.o \
apfelinit.o \
dd_test.o \

InitfION = \
rho.o \
ionpars.o \
rhonorm.o \
rhoxy.o \
rhoxycalc.o \
tpcalc.o \
tp.o \
rhoxyint.o \
betaion.o \
opacpcalc.o \
opacp.o \
opacpint.o \
opacpb.o \
opacpbcalc.o \
opacpbp.o \
opacpbpcalc.o \
opacpbint.o \
opacpbpint.o \
screencalc.o \
screen.o \
string.o \
s2qcdion.o \
s2qcdionp.o \
ioninit.o \
tpqcd.o \
tpqcdcalc.o \
tpqcdint.o \
gdrin.o \
gdr_int.o \
sudgam.o \
screen_01.o \

InitfUSER = \
alphasuser.o \
init.o \
initsud.o \
nf.o \
string.o \
hg.o \
inithg.o \
initpars.o \
calcop.o \
calcscreen.o \
opacityint.o \
screeningint.o \
screening.o \
opacity.o \
PDF.o \
PDFuser.o \
Sudakov.o \
inpdfuser.o \
dd.o \
Elastic.o \
sdcoh.o \
rann.o \
readscreen.o \
sdincoh.o \
F1F2.o \
F2ev.o \
rfit.o \
gd-fit-11.o \
CLAS.o \
apfelinit.o \

Helas = \
FFV2_0.o \
FFV2_1.o \
FFV2_2.o \
FFV2_3.o \
FFV1_0.o \
FFV1_1.o \
FFV1_2.o \
FFV1P0_3.o \
FFV3_0.o \
FFV3_1.o \
FFV3_2.o \
FFV3_3.o \
VVS1_0.o \
VVS1_3.o \
VVS1_3.o \
VVVV5_0.o \
FFV5_0.o \
FFV5_1.o \
FFV5_2.o \
FFV5_3.o \
VVVV2_0.o \
VVV1_0.o \
VVV1_1.o \
VVV1P0_1.o \
VVV1_2.o \
VVV1_3.o \
aloha_functions.o \
VVVV5_0.o \
VVVV2P0_1.o \
VVVV5_4.o \

Model = \
rw_para.o \
lha_read.o \
couplings.o \
couplings1.o \
couplings2.o \

MG = \
wwpol_match.o \
is_born.o \
matrix_uc.o \
matrix_ds.o \
matrix_ubcb.o \
matrix_dbsb.o \
matrix_dbub.o \
matrix_udb.o \
matrix_cub.o \
matrix_sdb.o \
matrix_dub.o \
matrix_du.o \
mgcross.o \
cs_SF_SD.o \
matrix_au.o \
matrix_aub.o \
matrix_ad.o \
matrix_adb.o \

#

sCODELHAi = $(Mainf) $(Mesf) $(EW) $(PDFsfLHA) $(Intf) $(Phasef) $(Subampsf) $(Survf) $(Userf) $(sPDFsf) $(Unwf) $(Varf) $(LbyLf) $(InitfION) $(Dissf) $(InitfLHA) $(Helas) $(Model) $(MG)

sCODELHA = $(patsubst %,$(OBJ_PATH)%,$(sCODELHAi))

iCODELHAi = $(InitfLHA) 
iCODELHA = $(patsubst %,$(OBJ_PATH)%,$(iCODELHAi))

iCODEIONi = $(InitfION) 
iCODEION = $(patsubst %,$(OBJ_PATH)%,$(iCODEIONi))

###########

all : init superchic superchicLib

$(sCODELHA): src/inc/head.f

src/inc/head.f: src/inc/head.f.inc
		cp src/inc/head.f.in src/inc/head.f
		gsed -i 's/@PROJECT_VERSION@/5.3/g' src/inc/head.f
		gsed -i 's/@RELEASE_DATE@/30.06.2024/g' src/inc/head.f

superchicLib: $(sCODELHA)
	$(FC) -L$(LHAPDFLIB) $(LIBFLAGlha) -mcmodel=large -shared -fPIC -o lib/libsuperchic.so $^
	ar rc lib/libsuperchic.a obj/*.o

superchic.o:	superchic.f
	$(FC) $(FFLAGS) -I$(INCPATH) -c  $< -o $@

init.o:	init.f
	$(FC) $(FFLAGS) -I$(INCPATH) -c  $< -o $@

$(OBJ_PATH)Elastic.o: Elastic.f
	$(FC) $(FFLAGS) -cpp -traditional -ffree-form -DDATA_PREFIX=`$$PREFIX/share/SuperChic` -I$(INCPATH) -c  $< -o $@

$(OBJ_PATH)gdrin.o: gdrin.f
	$(FC) $(FFLAGS) -cpp -traditional -ffree-form -DDATA_PREFIX=`$$PREFIX/share/SuperChic` -I$(INCPATH) -c  $< -o $@

$(OBJ_PATH)%.o: %.f
	$(FC) $(FFLAGS) -I$(INCPATH) -c  $< -o $@

superchic : $(OBJ_PATH)superchic.o $(sCODELHA)
	$(FC) $^ -L$(LHAPDFLIB) $(LIBFLAGlha) -o bin/$@

init : $(OBJ_PATH)init.o $(iCODELHA)
	$(FC) $^ -L$(LHAPDFLIB) $(LIBFLAGlha) -o bin/$@

clean:
	rm -f bin/init bin/superchic lib/lib* *.o $(OBJ_PATH)*.o

install: all
	@echo "Installing to $$PREFIX..."
	@mkdir -p $$PREFIX $$PREFIX/share/SuperChic
	@cp -r obj $$PREFIX/
	@cp -r bin $$PREFIX/
	@cp -r doc $$PREFIX/
	@cp -r lib $$PREFIX/
	@cp -r src $$PREFIX/
	@cp -r Cards $$PREFIX/
	@cp data/SplinesWithVariableKnots.dat  $$PREFIX/share/SuperChic
	@cp data/Carlos_106_440.dat  $$PREFIX/share/SuperChic
	@cp data/gampgamn.dat  $$PREFIX/share/SuperChic
	@cp data/Veyssiere_singleneut.dat  $$PREFIX/share/SuperChic
	@cp data/Muccifora.dat $$PREFIX/share/SuperChic
	@cp data/Lepretre25_103.dat  $$PREFIX/share/SuperChic
	@cp data/Caldwell.dat  $$PREFIX/share/SuperChic
	@cp data/Lepretre_25_103.dat  $$PREFIX/share/SuperChic
	@echo "Installation complete."
