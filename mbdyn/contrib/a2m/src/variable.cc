//                                 VARIABLE.CC                                  


#include <variable.h>
#include <output.h>


s_variable::s_variable(void) : Ic(0), Expression (new char[256]),
                               Value(0),
                               _Ic(N),_Expression(N),_Value(N) { }

s_variable::~s_variable(void) {}

inline const char* const s_variable::Gettype(void) const
{
   return "VARIABLE";
}

Boolean s_variable::Test(void)
{
   const int err_before=nerr;
   /* here is the test code */
   if (err_before != nerr) return Y; else return N;
}

ostream& s_variable::Print (ostream& out) const
{
   out << endl;
   out << "VARIABLE: " << label << endl;
   out << "      " << "Ic [" << _Ic << "] = " << Ic << endl
       << "      " << "Expression [" << _Expression << "] = " << Expression
       << endl
       << "      " << "Value [" << _Value << "] = " << Value
       << endl;
   out << endl;
   return out;
}

void s_variable::Translate (ostream& out)
{
   /* Here is the translation code */
   return;
}
