abundance_file=""
panman=""
output=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --abundance)
      abundance_file=$2
      shift 2
      ;;
    --panman)
      panman=$2
      shift 2
      ;;
    --output)
      output=$2
      shift 2
      ;;
    *)
      shift
      ;;
  esac
done

output_dir=$(dirname "$output")
output_base=$(basename "$output")
mkdir -p $output_dir
tmp_dir=$(mktemp -d "${output_dir}/${output_base}_tmp.XXXXXX")
trap 'rm -rf "$tmp_dir"' EXIT INT TERM HUP

fastas=()
for node in $(cut -f1 $abundance_file | tr ',' '\n' | grep '\S'); do
  clean_out=$(echo $node | tr '/' '_' | tr '|' '_').fa
  panmap $panman --dump-sequence $node -o $tmp_dir/$clean_out
  fastas+=($tmp_dir/$clean_out)
done

cat ${fastas[@]} > $tmp_dir/combined.fa

pangolin $tmp_dir/combined.fa --outfile $tmp_dir/combined.pangolin.csv

SCRIPT_DIR="$(dirname "${BASH_SOURCE[0]}")"
python3 $SCRIPT_DIR/assign_lineage.py $abundance_file $tmp_dir/combined.pangolin.csv > $output