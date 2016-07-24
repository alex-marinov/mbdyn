#include "mbconfig.h"
#include "solver.h"
#include "invsolver.h"

struct mb_sol_wrap_t {
	Table *pT;
	MathParser *pMP;
	std::ifstream streamIn;
	InputStream *pIn;
	MBDynParser *pHP;
	Solver *pS;
};

static mb_sol_wrap_t *
mb_sol_create_(const char *sIn, const char *sOut)
{
	mb_sol_wrap_t *pSW = new mb_sol_wrap_t;
	pSW->pT = new Table(true);
	pSW->pMP = new MathParser(*pSW->pT, false);
	pSW->streamIn.open(sIn);
	pSW->pIn = new InputStream(pSW->streamIn);

	return pSW;
}

extern "C" void *
mb_sol_create(const char *sIn, const char *sOut)
{
	mb_sol_wrap_t *pSW = mb_sol_create_(sIn, sOut);
	pSW->pS = new Solver(*pSW->pHP, sIn, sOut, false);
	return (void *)pSW;
}

extern "C" void *
mb_sol_create_inv(const char *sIn, const char *sOut)
{
	mb_sol_wrap_t *pSW = mb_sol_create_(sIn, sOut);
	pSW->pS = new InverseSolver(*pSW->pHP, sIn, sOut, false);
	return (void *)pSW;
}

extern "C" int
mb_sol_prepare(void *p)
{
	mb_sol_wrap_t *pSW = (mb_sol_wrap_t *)p;
	if (!pSW->pS->Prepare()) {
		return -1;
	}

	return 0;
}

extern "C" int
mb_sol_start(void *p)
{
	mb_sol_wrap_t *pSW = (mb_sol_wrap_t *)p;
	if (!pSW->pS->Start()) {
		return -1;
	}

	return 0;
}

extern "C" int
mb_sol_advance(void *p)
{
	mb_sol_wrap_t *pSW = (mb_sol_wrap_t *)p;
	if (!pSW->pS->Advance()) {
		return -1;
	}

	return 0;
}

extern "C" int
mb_sol_destroy(void *p)
{
	mb_sol_wrap_t *pSW = (mb_sol_wrap_t *)p;
	delete pSW->pS;
	delete pSW->pIn;
	delete pSW->pMP;
	delete pSW->pT;
	delete pSW;

	return 0;
}

extern "C" int
mb_sol_setbufin(void *p, unsigned uLabel, integer iSize, doublereal *pdBuf)
{
	mb_sol_wrap_t *pSW = (mb_sol_wrap_t *)p;
	DataManager *pDM = pSW->pS->pGetDataManager();
	pDM->SetBufInRaw(uLabel, iSize, pdBuf);

	return 0;
}

extern "C" int
mb_sol_setbufout(void *p, unsigned uLabel, integer iSize, doublereal *pdBuf)
{
	mb_sol_wrap_t *pSW = (mb_sol_wrap_t *)p;
	DataManager *pDM = pSW->pS->pGetDataManager();
	pDM->SetBufOutRaw(uLabel, iSize, pdBuf);

	return 0;
}
