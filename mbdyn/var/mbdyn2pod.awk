# Skip undesired labels
function skip_label(l)
{
	return 0;
}

BEGIN {
	FirstStep = 1;
	FirstLabel = -1;
} 
{
	# if first label, end line
	if ($1 == FirstLabel) {
		printf("\n");
	}

	# set first label
	if (FirstStep) {
		FirstLabel = $1
		FirstStep = 0;
	}

	# skip undesired labels (modify the body of the function at will)
	if (skip_label($1)) {
		next;
	}

	# collect values
	if (last == 0) {
		l = NF;
	} else {
		l = last;
	}

	for (j = 2; j <= l; j++) {
		printf(" %13.6e", $(j));
	}
}
END {
	# end last line
	printf("\n");
}
