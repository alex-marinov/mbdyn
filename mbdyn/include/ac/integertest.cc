#include "mbconfig.h"

#include "ac/f2c.h"

template<class T> class a{
public:
	static int x;
};

#ifdef INTEGER_TEST
int a<int>::x = 0;
int a<integer>::x = 0;
#endif /* INTEGER_TEST */

int
main(void)
{
#ifdef INTEGER_TEST
	a<int> xint;
	a<integer> xinteger;

	if (&xint.x == &xinteger.x) {
		return 0;
	}
#endif /* INTEGER_TEST */

	return 1;
}

