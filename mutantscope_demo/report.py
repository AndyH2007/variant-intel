from reportlab.lib.pagesizes import A4
from reportlab.pdfgen import canvas

def write_report(path, title="MutantScope Report", lines=None):
    c = canvas.Canvas(path, pagesize=A4)
    w, h = A4
    y = h - 72
    c.setFont("Helvetica-Bold", 16)
    c.drawString(72, y, title)
    y -= 24
    c.setFont("Helvetica", 11)
    if lines:
        for line in lines:
            c.drawString(72, y, str(line))
            y -= 14
            if y < 72:
                c.showPage(); y = h - 72
    c.showPage()
    c.save()
