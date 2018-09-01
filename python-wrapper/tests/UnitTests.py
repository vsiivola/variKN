#!/usr/bin/env python3

import unittest, sys, os
sys.path.append("../../lib/python")
import varikn

class TestVarigramTraining(unittest.TestCase):

    def setUp(self):
        self.data = "../../tests/data/abx20.txt"
        self.vocab = "../../tests/data/vocab.txt"
        self.lm = "../../tests/data/a2.arpa"

    def test_training(self):
        vg = varikn.VarigramTrainer(False, False)
        vg.set_datacost_scale(1.0);
        vg.set_datacost_scale2(2.0);
        vg.initialize(self.data, 100, 0, 9999999,
                      self.data, "", False, "");
        vg.grow(1);
        vg.write_file(self.lm, True);
        assert(os.path.isfile(self.lm))

class TestInterpolation(unittest.TestCase):

    def setUp(self):
        self.lm1 = "../../tests/data/a.arpa"
        self.lm2 = "../../tests/data/b.arpa"
        self.data = "../../tests/data/abx20.txt"

    def test_interpolation(self):
        tg = varikn.InterTreeGram([self.lm1, self.lm2], [0.6, 0.4])
        lm = varikn.Perplexity(tg, "", "", "", "")
        lm.set_init_hist(0)
        lpsum = 0.0
        with open(self.data) as fobj:
            for line in fobj:
                for token in line.split():
                    lpsum += lm.token_logprob(token)
                lm.clear_history()
        print(lpsum)

if __name__ == '__main__':
    unittest.main()
