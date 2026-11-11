import argparse
import gzip
import os
import re
import shutil
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.mast import MastMissions, Observations


def sanitize_folder_name(name):
  clean_name = re.sub(r"[^\w\-]", "_", name.strip())
  return re.sub(r"_+", "_", clean_name).strip("_")


def get_dataset_ids(obs_table):
  for col in [
      "iue_data_id",
      "data_set_name",
      "sci_data_set_name",
      "obs_id",
      "dataset_id",
  ]:
    if col in obs_table.colnames:
      return [str(val).strip() for val in obs_table[col]]
  raise KeyError(
      f"Could not find dataset ID column. Available columns: {obs_table.colnames}"
  )


def filter_mxlo_gz_table(data_products):
  file_col = None
  for col in ["filename", "productFilename", "dataURI", "url"]:
    if col in data_products.colnames:
      file_col = col
      break

  if file_col is None:
    print(f"Available product columns: {data_products.colnames}")
    return data_products[0:0]

  # Match only files ending with .mxlo.gz (e.g., lwr14197.mxlo.gz)
  mask = [
      str(fn).strip().lower().endswith(".mxlo.gz")
      for fn in data_products[file_col]
  ]
  return data_products[mask]


def extract_and_flatten_downloads(download_dir):
  extracted_files = []
  for root, dirs, files in os.walk(download_dir, topdown=False):
    for filename in files:
      if filename.lower().endswith(".mxlo.gz"):
        src_gz_path = os.path.join(root, filename)
        out_filename = filename[:-3]
        dst_mxlo_path = os.path.join(download_dir, out_filename)

        # Decompress directly into download_dir
        with gzip.open(src_gz_path, "rb") as f_in:
          with open(dst_mxlo_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

        # Remove the original .mxlo.gz file
        os.remove(src_gz_path)
        extracted_files.append(dst_mxlo_path)

    # Remove empty subdirectories left behind by astroquery
    if os.path.abspath(root) != os.path.abspath(download_dir):
      try:
        os.rmdir(root)
      except OSError:
        pass

  return extracted_files


def download_iue_mxlo(
    star_name,
    search_radius_arcsec=180.0,
    download_dir=None,
    manual_coords=None,
    resolver="SESAME",
):
  if download_dir is None:
    download_dir = os.path.join(".", sanitize_folder_name(star_name))

  iue_missions = MastMissions(mission="iue")

  # 1. Coordinate Lookup
  if manual_coords:
    print(f"Using manual coordinates: {manual_coords}")
    coords = SkyCoord(manual_coords, unit=(u.hourangle, u.deg))
  else:
    print(f"Resolving coordinates for '{star_name}' using {resolver}...")
    try:
      if resolver.upper() in ("SESAME", "SESEME"):
        coords = SkyCoord.from_name(star_name)
      else:
        coords = iue_missions.resolve_object(star_name, resolver=resolver)
    except Exception as e:
      print(f"Error resolving '{star_name}': {e}")
      return

  print(
      f"Coordinates used: RA = {coords.ra.deg:.5f} deg, Dec ="
      f" {coords.dec.deg:.5f} deg"
  )

  # 2. Query IUE mission metadata near the coordinates
  print(f"\nSearching MAST IUE mission within {search_radius_arcsec}\"...")
  obs_table = iue_missions.query_region(
      coordinates=coords, radius=search_radius_arcsec * u.arcsec
  )

  if len(obs_table) == 0:
    print("No IUE observations found for this target.")
    return

  print(f"Found {len(obs_table)} IUE observation(s).")
  dataset_ids = get_dataset_ids(obs_table)
  print(f"Dataset IDs found: {dataset_ids}")

  # 3. Retrieve and filter .mxlo.gz products
  try:
    data_products = iue_missions.get_unique_product_list(dataset_ids)
    print(f"Total products found: {len(data_products)}")

    filtered_products = filter_mxlo_gz_table(data_products)
    print(f"Filtered down to {len(filtered_products)} .mxlo.gz file(s).")

    if len(filtered_products) > 0:
      manifest = iue_missions.download_products(
          filtered_products, download_dir=download_dir
      )
      extracted = extract_and_flatten_downloads(download_dir)
      print(
          f"\nDownload and extraction complete. {len(extracted)} .mxlo file(s)"
          f" saved directly in: {download_dir}"
      )
      return manifest
    else:
      print("\nNo .mxlo.gz files matched. All available products:")
      data_products.pprint(max_width=-1)

  except Exception as e:
    print(f"MastMissions product lookup failed ({e}). Falling back to CAOM...")
    caom_ids = [did.lower() for did in dataset_ids]
    caom_obs = Observations.query_criteria(
        obs_collection="IUE", obs_id=caom_ids
    )
    data_products = Observations.get_product_list(caom_obs)
    filtered_products = filter_mxlo_gz_table(data_products)
    print(f"Filtered down to {len(filtered_products)} .mxlo.gz file(s).")

    if len(filtered_products) > 0:
      manifest = Observations.download_products(
          filtered_products, download_dir=download_dir
      )
      extracted = extract_and_flatten_downloads(download_dir)
      print(
          f"\nDownload and extraction complete. {len(extracted)} .mxlo file(s)"
          f" saved directly in: {download_dir}"
      )
      return manifest
    else:
      print("No .mxlo.gz files found for this target.")


if __name__ == "__main__":
  parser = argparse.ArgumentParser(
      description=(
          "Download and extract low-dispersion (.mxlo.gz) IUE spectra using"
          " MastMissions."
      )
  )
  parser.add_argument(
      "target",
      type=str,
      help="Name or catalog ID of the target star (e.g., 'sk-66 19')",
  )
  parser.add_argument(
      "-r",
      "--radius",
      type=float,
      default=60.0,
      help="Search radius in arcseconds (default: 60.0)",
  )
  parser.add_argument(
      "-o",
      "--outdir",
      type=str,
      default=None,
      help=(
          "Custom output directory (defaults to a folder named after the star)"
      ),
  )
  parser.add_argument(
      "--resolver",
      type=str,
      default="SESAME",
      choices=["SESAME", "SIMBAD", "NED"],
      help="Name resolver to use (default: SESAME)",
  )
  parser.add_argument(
      "-c",
      "--coords",
      type=str,
      default=None,
      help=(
          "Manual RA/Dec coordinates to bypass resolver (e.g., '18h36m56.3s"
          " +38d47m01s')"
      ),
  )

  args = parser.parse_args()
  download_iue_mxlo(
      args.target,
      search_radius_arcsec=args.radius,
      download_dir=args.outdir,
      manual_coords=args.coords,
      resolver=args.resolver,
  )
