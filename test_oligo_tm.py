import csv
import json

from oligo_tm import calculate_metrics, main, process_batch


def test_oligo_tm_single(capsys):
    result = calculate_metrics("ATCGATCGATCGATCGATCG", na_mM=50, mg_mM=1.5)
    assert result["tm_mg_corrected_C"] == 59.72
    assert main(["single", "--sequence", "ATCGATCGATCGATCGATCG"]) == 0
    assert json.loads(capsys.readouterr().out)["sequence"] == "ATCGATCGATCGATCGATCG"


def test_oligo_tm_bad_sequence(capsys):
    assert main(["single", "--sequence", "NOT-DNA"]) == 2
    assert "error" in json.loads(capsys.readouterr().out)


def test_oligo_tm_batch(tmp_path):
    source = tmp_path / "in.csv"
    target = tmp_path / "out.csv"
    source.write_text("id,sequence\nPrimer_1,ATCGATCGATCGATCGATCG\nPrimer_2,INVALID\n", encoding="utf-8")
    assert process_batch(str(source), str(target)) == 2
    with target.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    assert rows[0]["tm_mg_corrected_C"] == "59.72"
    assert rows[1]["error"]
    assert rows[1]["id"] == "Primer_2"
