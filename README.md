VariKN language modeling toolkit provides tools for training n-gram
language models. Amongst the supported methods are:

- Absolute discounting
- Kneser-Ney smoothing
- Revised Kneser pruning
- Kneser-Ney growing

The descriptions of Revised Kneser pruning and Kneser-Ney growing can
be found in the paper Vesa Siivola, Teemu Hirsimäki and Sami Virpioja,
"On Growing and Pruning Kneser-Ney Smoothed N-Gram Models", IEEE
Transactions on Speech, Audio and Language Processing,
15(5):1617-1624, 2007.

The package provides a accurate pruning for Kneser-Ney smoothed
models. Also, it is possible to train a very high-order n-gram models
with the growing algorithm. The models can be output to arpa lm
format, which is compatible with most common other tools in the
field. You can look at the paper Vesa Siivola, Mathias Creutz and
Mikko Kurimo: "Morfessor and VariKN machine learning tools for speech
and language technology", Proceedings of the 8th International
Conference on Speech Communication and Technology (INTERSPEECH'07),
2007 for guidelines on typical use.

For help on installing the toolkit, read [INSTALL.md](INSTALL.md).

For help on using the built commands, read [commands.html](commands.html).

If it looks like there is interest in using the toolkit but the
documentation is inadequate, let me know the main sticking points and
I'll try to do something about it. If you are interested in a well
supported toolkit that provides a wide range of functionality, you
might want to look at the [srilm toolkit](http://www.speech.sri.com/projects/srilm/).
