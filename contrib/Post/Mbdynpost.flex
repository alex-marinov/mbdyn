%option c++
%option warn
%option prefix="Mbdynpost"
%option noyywrap
%option debug

%{

#include "Post.hh"

%}


ws              [ \t]*

output_frequency_tok		"output frequency:"
step_tok			"Step"

alpha		[A-Za-z]
dig		[0-9]
comma      	","
duepunti  	":"
dueduepunti  	"::"
name		({alpha}|{dig})+
fileextension	{alpha}+
label		{dig}+

comment		"#"


%%

{output_frequency_tok}			{return OUTPUT_FREQUENCY_TOK;}
{step_tok}				{return STEP_TOK;}
{comment}				{	/* skip c++-style comments */
						int c;
						while((c = yyinput()) != 0) {
							if(c == '\n') {
								break;
							}
						}
					}
{fileextension}				{return FILE_EXTENSION_TOK;}
{label}					{return LABEL_TOK;}
{comma}					{return COMMA_TOK;}
{duepunti}				{return DUEPUNTI_TOK;}
{dueduepunti}				{return DUEDUEPUNTI_TOK;}
\n					{}
.					{}
<<EOF>>					{return EOF_TOK;}
