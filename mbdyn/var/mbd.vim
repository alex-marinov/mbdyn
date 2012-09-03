" $Header$
" Vim syntax file
" Language:	MBDyn
" Maintainer:	Pierangelo Masarati <masarati@aero.polimi.it>
" Last Change:	2008 Feb 29
" based on c syntax by Bram Moolenaar <Bram@vim.org>

" For version 5.x: Clear all syntax items
" For version 6.x: Quit when a syntax file was already loaded
if version < 600
  syntax clear
elseif exists("b:current_syntax")
  finish
endif

syn keyword	mbdTodo	contained TODO FIXME XXX

" cCommentGroup allows adding matches for special things in comments
syn cluster	mbdCommentGroup	contains=mbdTodo

" String and Character constants
" Highlight special characters (those which have a backslash) differently
syn match	mbdSpecial	display contained "\\\(x\x\+\|\o\{1,3}\|.\|$\)"
syn region	mbdString	start=+L\="+ skip=+\\\\\|\\"+ end=+"+ contains=mbdSpecial,@Spell

syn match	mbdCharacter	"L\='[^\\]'"
syn match	mbdCharacter	"L'[^']*'" contains=mbdSpecial

"integer number, or floating point number without a dot and with "f".
syn case ignore
" A bunch of useful MBDyn keywords
"syn region mbdBlock start=/begin\s*:/ms=e+1 end=/;/me=s-1 contains=@mbdBlockStmt
"syn region mbdBlock start=/end\s*:/ms=e+1 end=/;/me=s-1 contains=@mbdBlockStmt
"syn region mbdBlockStmt contained start=/:/ms=e+1 end=/;/me=s-1
syn region mbdBlock start=/begin\s*:/ end=/;/me=s-1 contains=mbdBlockStmt2
syn region mbdBlock start=/end\s*:/ end=/;/me=s-1 contains=mbdBlockStmt2
syn region mbdBlockStmt2 start=/:/ end=/;/me=s-1 contains=mbdBlockStmt contained 
syn region mbdBlockStmt start=/:/ms=e+1 end=/;/me=s-1 contained 

" MBDyn description
syn match	mbdDescription	"[_a-zA-Z][_a-zA-Z0-9[:space:]]*:"me=e-1 contains=ALLBUT,@mbdBlock
syn match	mbdDescription	";[[:space:]\n]*[_a-zA-Z][_a-zA-Z0-9[:space:]]*;"ms=s+1,me=e-1

syn match	mbdNumbers	display transparent "\<\d\|\.\d" contains=mbdNumber,mbdFloat,mbdOctalError,mbdOctal
" Same, but without octal error (for comments)
syn match	mbdNumbersCom	display contained transparent "\<\d\|\.\d" contains=mbdNumber,mbdFloat,mbdOctal
syn match	mbdNumber	display contained "\d\+\(u\=l\{0,2}\|ll\=u\)\>"
"hex number
syn match	mbdNumber	display contained "0x\x\+\(u\=l\{0,2}\|ll\=u\)\>"
" Flag the first zero of an octal number as something special
syn match	mbdOctal	display contained "0\o\+\(u\=l\{0,2}\|ll\=u\)\>" contains=mbdOctalZero
syn match	mbdOctalZero	display contained "\<0"
syn match	mbdFloat	display contained "\d\+f"
"floating point number, with dot, optional exponent
syn match	mbdFloat	display contained "\d\+\.\d*\(e[-+]\=\d\+\)\=[fl]\="
"floating point number, starting with a dot, optional exponent
syn match	mbdFloat	display contained "\.\d\+\(e[-+]\=\d\+\)\=[fl]\=\>"
"floating point number, without dot, with exponent
syn match	mbdFloat	display contained "\d\+e[-+]\=\d\+[fl]\=\>"

" flag an octal number with wrong digits
syn match	mbdOctalError	display contained "0\o*[89]\d*"
syn case match

if exists("c_comment_strings")
  " A comment can contain cString, cCharacter and cNumber.
  " But a "*/" inside a cString in a cComment DOES end the comment!  So we
  " need to use a special type of cString: cCommentString, which also ends on
  " "*/", and sees a "*" at the start of the line as comment again.
  " Unfortunately this doesn't very well work for // type of comments :-(
  syntax match	mbdCommentSkip	contained "^\s*\*\($\|\s\+\)"
  syntax region mbdCommentString	contained start=+L\=\\\@<!"+ skip=+\\\\\|\\"+ end=+"+ end=+\*/+me=s-1 contains=mbdSpecial,mbdCommentSkip
  syntax region mbdComment2String	contained start=+L\=\\\@<!"+ skip=+\\\\\|\\"+ end=+"+ end="$" contains=mbdSpecial
  syntax region mbdCommentL	start="//" skip="\\$" end="$" keepend contains=@mbdCommentGroup,mbdComment2String,mbdCharacter,mbdNumbersCom,mbdSpaceError,@Spell
  syntax region mbdComment	matchgroup=mbdCommentStart start="/\*" end="\*/" contains=@mbdCommentGroup,mbdCommentStartError,mbdCommentString,mbdCharacter,mbdNumbersCom,mbdSpaceError,@Spell
else
  syn region	mbdCommentL	start="//" skip="\\$" end="$" keepend contains=@mbdCommentGroup,mbdSpaceError,@Spell
  syn region	mbdComment	matchgroup=mbdCommentStart start="/\*" end="\*/" contains=@mbdCommentGroup,mbdCommentStartError,mbdSpaceError,@Spell
endif
syntax match	mbdCommentError	display "\*/"
syntax match	mbdCommentStartError display "/\*"me=e-1 contained

" Comments
"=========
syn cluster    mbdShCommentGroup	contains=mbdShTodo,@Spell
syn keyword    mbdShTodo	contained	TODO
syn match      mbdShComment	"#.*$" contains=@mbdShCommentGroup

syn keyword	mbdType	bool integer real string const ifndef

syn keyword	mbdConstant pi

" Define the default highlighting.
" For version 5.7 and earlier: only when not done already
" For version 5.8 and later: only when an item doesn't have highlighting yet
if version >= 508 || !exists("did_mbd_syn_inits")
  if version < 508
    let did_mbd_syn_inits = 1
    command -nargs=+ HiLink hi link <args>
  else
    command -nargs=+ HiLink hi def link <args>
  endif

  HiLink mbdFormat		mbdSpecial
  HiLink mbdCppString		mbdString
  HiLink mbdCommentL		mbdComment
  HiLink mbdCommentStart	mbdComment
  HiLink mbdLabel		Label
  HiLink mbdUserLabel		Label
  HiLink mbdConditional		Conditional
  HiLink mbdDescription		Conditional
  HiLink mbdRepeat		Repeat
  HiLink mbdCharacter		Character
  HiLink mbdSpecialCharacter	mbdSpecial
  HiLink mbdNumber		Number
  HiLink mbdOctal		Number
  HiLink mbdOctalZero		PreProc	 " link this to Error if you want
  HiLink mbdFloat		Float
  HiLink mbdOctalError		mbdError
  HiLink mbdParenError		mbdError
  HiLink mbdErrInParen		mbdError
  HiLink mbdErrInBracket	mbdError
  HiLink mbdCommentError	mbdError
  HiLink mbdCommentStartError	mbdError
  HiLink mbdSpaceError		mbdError
  HiLink mbdSpecialError	mbdError
  HiLink mbdOperator		Operator
  HiLink mbdInclude		Include
  HiLink mbdPreProc		PreProc
  HiLink mbdDefine		Macro
  HiLink mbdIncluded		mbdString
  HiLink mbdError		Error
  HiLink mbdStatement		Statement
  HiLink mbdBlock		PreProc
  HiLink mbdBlockStmt		Identifier
  HiLink mbdType		Type
  HiLink mbdConstant		Constant
  HiLink mbdCommentString	mbdString
  HiLink mbdComment2String	mbdString
  HiLink mbdCommentSkip		mbdComment
  HiLink mbdString		String
  HiLink mbdComment		Comment
  HiLink mbdShComment		Comment
  HiLink mbdSpecial		SpecialChar
  HiLink mbdTodo		Todo

  delcommand HiLink
endif

let b:current_syntax = "mbd"

" vim: ts=8
