"""Run without protein kernels: python -m unittest discover -s . -p test_receiver_ranking.py."""
import unittest
import json
import tempfile
from pathlib import Path
import numpy as np
from ranking import validated_indices, receiver_rms


class ReceiverRankingTests(unittest.TestCase):
    def test_nonconsecutive_index_zero_is_a_receiver(self):
        indices = np.array([0, 3, 6], dtype=np.int64)
        mask = np.array([True, False, False, True, False, False, True])
        np.testing.assert_array_equal(validated_indices(indices, 7, kind='indices'), indices)
        np.testing.assert_array_equal(validated_indices(mask, 7, kind='mask'), indices)
        self.assertFalse(np.array_equal(np.flatnonzero(indices), indices))

    def test_old_conversion_changes_ranking_and_signal_floor(self):
        response = np.zeros((8, 8))
        response[:6, [1]] = np.arange(6, 0, -1).reshape(-1, 1) * 0.1
        response[:6, [3]] = np.arange(1, 7).reshape(-1, 1) * 0.0001
        response[:6, [7]] = np.arange(1, 7).reshape(-1, 1) * 0.0001
        receiver = np.array([3, 7])
        corrected = receiver_rms(response, receiver)
        wrong = receiver_rms(response, np.flatnonzero(receiver))
        expected = np.array([sum(response[i, j] ** 2 for j in [3, 7]) / 2 for i in range(8)]) ** 0.5
        np.testing.assert_allclose(corrected, expected, rtol=0, atol=1e-16)
        np.testing.assert_array_equal(np.argsort(-corrected)[:5], [5, 4, 3, 2, 1])
        np.testing.assert_array_equal(np.argsort(-wrong)[:5], [0, 1, 2, 3, 4])
        self.assertLessEqual(corrected.max(), .002)
        self.assertGreater(wrong.max(), .002)

    def test_boolean_and_integer_encodings_are_not_interchangeable(self):
        with self.assertRaises(TypeError):
            validated_indices(np.array([False, True]), 2, kind='indices')
        with self.assertRaises(TypeError):
            validated_indices(np.array([0, 1]), 2, kind='mask')
        with self.assertRaises(TypeError):
            receiver_rms(np.ones((2, 2)), np.array([False, True]))

    def test_all_four_saved_protein_receiver_scores(self):
        root = Path(__file__).resolve().parent
        for target in ['kras', 'abl', 'myc', 'myh7']:
            with self.subTest(target=target), np.load(root/'fixtures'/f'{target}-tau.npz') as fixture:
                ids = fixture['canonical']
                candidates = validated_indices(fixture['candidate'], len(ids), kind='mask')
                receivers = validated_indices(fixture['receiver'], len(ids), kind='indices')
                row = json.loads((root/'results'/target/'analysis.json').read_text())['results']['selected_times']['rows'][1]
                for key, field in [('reference', 'reference_top5'), ('reduced', 'reduced_top5')]:
                    scores = receiver_rms(fixture[key], receivers)
                    ordered = candidates[np.lexsort((ids[candidates], -scores[candidates]))]
                    self.assertEqual(ids[ordered[:5]].tolist(), row[field])
                    wrong_scores = receiver_rms(fixture[key], np.flatnonzero(fixture['receiver']))
                    wrong_order = candidates[np.lexsort((ids[candidates], -wrong_scores[candidates]))]
                    self.assertNotEqual(ids[wrong_order[:5]].tolist(), row[field])
                self.assertGreater(scores[candidates].max(), .002)

    def test_api_refuses_existing_analysis_before_reading_inputs(self):
        from analyze import run
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            analysis = output/'results/kras/analysis.json'
            analysis.parent.mkdir(parents=True)
            analysis.write_text('unchanged')
            with self.assertRaises(FileExistsError):
                run('kras', output/'missing-inputs', output/'missing-reference', output)
            self.assertEqual(analysis.read_text(), 'unchanged')

    def test_invalid_selectors_fail_explicitly(self):
        for values in [np.array([-1]), np.array([3]), np.array([1, 1]), np.array([], dtype=int), np.array([[1]])]:
            with self.subTest(values=values), self.assertRaises(ValueError):
                validated_indices(values, 3, kind='indices')
        for values in [np.array([1.0]), np.array(['1']), np.array([1], dtype=object)]:
            with self.subTest(values=values), self.assertRaises(TypeError):
                validated_indices(values, 3, kind='indices')
        with self.assertRaises(ValueError):
            validated_indices(np.array([True]), 3, kind='mask')


if __name__ == '__main__':
    unittest.main()
