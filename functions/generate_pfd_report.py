from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
from reportlab.lib.units import inch
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Image


def create_pdf_report(df_list, figs, filename, text_elements):
    """
    Creates a PDF report with the given data and saves it to a file.

    Args:
        df_list (list): A list of DataFrame objects containing the data to be included in the report.
        figs (list): A list of figure objects to be included in the report.
        filename (str): The name of the file to save the PDF report to.
        text_elements (dict): A dictionary containing various text elements to be included in the report.

    Returns:
        None
    """
    # Create a new PDF document
    doc = SimpleDocTemplate(filename, pagesize=letter)

    # Create a list of flowables to add to the document
    flowables = []

    # Add a title to the document
    title_style = getSampleStyleSheet()["Title"]
    flowables.append(Paragraph("PDF Report " + filename + ".pdf", title_style))

    # Add an empty paragraph to the flowables list
    flowables.append(Paragraph("", text_style))

    # Add a paragraph of text to the document
    text_style = getSampleStyleSheet()["Normal"]
    flowables.append(
        Paragraph(
            "Report of results from UTOPIA",
            text_style,
        )
    )
    # Add an empty paragraph to the flowables list
    flowables.append(Paragraph("", text_style))

    # Add a new paragraph of text to the document
    new_text = (
        "Emission of plastic particles of "
        + str(text_elements["imput_flow_g_s"])
        + " g per second into the "
        + text_elements["recieving_compartment"]
    )
    flowables.append(Paragraph(new_text, text_style))
    flowables.append(
        Paragraph(
            "Particles from: " + text_elements["particle_emissions_form"], text_style
        )
    )
    flowables.append(
        Paragraph(
            "Partciles size: "
            + str(text_elements["particle_emissions_size_nm"])
            + " nm",
            text_style,
        )
    )
    flowables.append(
        Paragraph(
            "Particles density: "
            + str(text_elements["plastic_density_kg_m3"])
            + " kg-m3",
            text_style,
        )
    )

    # Add an empty paragraph to the flowables list
    flowables.append(Paragraph("", text_style))

    # Define a style for the subsection title
    subsection_style = getSampleStyleSheet()["Heading2"]

    # Add the subsection title to the document
    flowables.append(Paragraph("Results by compartment", subsection_style))

    # Add content for the subsection

    # Add figures to the document
    # for fig in figs:
    #     flowables.append(Image(fig, width=6 * inch))

    # Add a table to the document
    for df in df_list:
        table_data = [list(df.columns)] + df.values.tolist()
        table = Table(table_data)
        table.setStyle(
            TableStyle(
                [
                    ("BACKGROUND", (0, 0), (-1, 0), colors.grey),
                    ("TEXTCOLOR", (0, 0), (-1, 0), colors.whitesmoke),
                    (
                        "LINEBELOW",
                        (0, 0),
                        (-1, 0),
                        1,
                        colors.black,
                    ),  # Add a line below the header row
                    ("ALIGN", (0, 0), (-1, 0), "CENTER"),
                    ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                    ("FONTSIZE", (0, 0), (-1, 0), 14),
                    ("BOTTOMPADDING", (0, 0), (-1, 0), 12),
                    ("BACKGROUND", (0, 1), (-1, -1), colors.beige),
                    ("TEXTCOLOR", (0, 1), (-1, -1), colors.black),
                    ("ALIGN", (0, 1), (-1, -1), "CENTER"),
                    ("FONTNAME", (0, 1), (-1, -1), "Helvetica"),
                    ("FONTSIZE", (0, 1), (-1, -1), 12),
                    ("BOTTOMPADDING", (0, 1), (-1, -1), 6),
                    ("GRID", (0, 0), (-1, -1), 1, colors.black),
                ]
            )
        )
        flowables.append(table)

    # Build the PDF document
    doc.build(flowables)


# # Example usage
# import pandas as pd
# import matplotlib.pyplot as plt

# # Create a sample dataframe
# df = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6], "C": [7, 8, 9]})

# # Create a sample figure
# fig1, ax1 = plt.subplots()
# ax1.plot([1, 2, 3], [4, 5, 6])
# fig1.savefig("fig1.png")

# # Create another sample figure
# fig2, ax2 = plt.subplots()
# ax2.plot([1, 2, 3], [7, 8, 9])
# fig2.savefig("fig2.png")

# # Generate the PDF report
create_pdf_report(df_list, figs, filename, text_elements)


# def create_pdf_compartments_subsections(comp,
