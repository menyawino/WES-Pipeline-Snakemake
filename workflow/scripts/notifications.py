#!/usr/bin/env python3
"""
Pipeline Notification Module.
Handles email notifications upon workflow completion (onsuccess) or failure (onerror).
Supports SMTP authentication, STARTTLS, SSL, and local mailer fallbacks with rich HTML reports.
"""

import os
import sys
import smtplib
import logging
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from datetime import datetime

logger = logging.getLogger(__name__)

def build_html_report(status, pipeline_name, outdir, samples, log_file=None, error_details=None):
    """Generate a clean, aesthetic HTML email report body."""
    is_success = (status.lower() == "success")
    status_color = "#10b981" if is_success else "#ef4444"
    status_icon = "✓" if is_success else "✗"
    status_text = "Completed Successfully" if is_success else "Execution Failed"
    
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    sample_list_html = "".join(f"<li><code>{s}</code></li>" for s in samples) if samples else "<li>None specified</li>"
    
    acmg_dashboard = os.path.join(outdir, "analysis", "007_annotation", "cohort_acmg_dashboard.html")
    variant_dashboard = os.path.join(outdir, "analysis", "010_summary", "cohort_variant_dashboard.html")
    
    links_html = ""
    if is_success:
        links_html = f"""
        <div style="margin-top: 20px; padding: 15px; background: #f0fdf4; border-left: 4px solid #10b981; border-radius: 4px;">
            <h4 style="margin: 0 0 10px 0; color: #166534;">Generated Reports & Dashboards</h4>
            <ul style="margin: 0; padding-left: 20px; color: #14532d;">
                <li><strong>ACMG Dashboard:</strong> <code>{acmg_dashboard}</code></li>
                <li><strong>Cohort Variant Dashboard:</strong> <code>{variant_dashboard}</code></li>
            </ul>
        </div>
        """

    error_html = ""
    if not is_success and (error_details or log_file):
        error_html = f"""
        <div style="margin-top: 20px; padding: 15px; background: #fef2f2; border-left: 4px solid #ef4444; border-radius: 4px;">
            <h4 style="margin: 0 0 10px 0; color: #991b1b;">Error Details</h4>
            <p style="margin: 0; color: #7f1d1d;"><strong>Log File:</strong> <code>{log_file or 'N/A'}</code></p>
            {f'<pre style="background: #ffffff; padding: 10px; border-radius: 4px; overflow-x: auto; font-size: 12px; color: #b91c1c;">{error_details}</pre>' if error_details else ''}
        </div>
        """

    html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="utf-8">
        <style>
            body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif; line-height: 1.6; color: #1f2937; margin: 0; padding: 0; background-color: #f9fafb; }}
            .container {{ max-width: 650px; margin: 20px auto; background: #ffffff; border-radius: 8px; box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1); overflow: hidden; border: 1px solid #e5e7eb; }}
            .header {{ background: #111827; color: #ffffff; padding: 24px; text-align: center; }}
            .header h1 {{ margin: 0; font-size: 20px; font-weight: 700; letter-spacing: 0.5px; }}
            .header p {{ margin: 5px 0 0 0; color: #9ca3af; font-size: 13px; }}
            .badge {{ display: inline-block; padding: 6px 14px; border-radius: 9999px; font-weight: 600; font-size: 13px; color: #ffffff; background-color: {status_color}; margin-top: 15px; }}
            .content {{ padding: 24px; }}
            .meta-table {{ width: 100%; border-collapse: collapse; margin-top: 15px; }}
            .meta-table td {{ padding: 8px 12px; border-bottom: 1px solid #f3f4f6; font-size: 13px; }}
            .meta-table td.label {{ font-weight: 600; color: #4b5563; width: 30%; }}
            .meta-table td.value {{ color: #111827; }}
            .footer {{ background: #f3f4f6; padding: 15px; text-align: center; font-size: 12px; color: #6b7280; border-top: 1px solid #e5e7eb; }}
            code {{ background: #f3f4f6; padding: 2px 5px; border-radius: 3px; font-family: monospace; font-size: 12px; }}
        </style>
    </head>
    <body>
        <div class="container">
            <div class="header">
                <h1>{pipeline_name}</h1>
                <p>Automated Execution Notification</p>
                <div class="badge">{status_icon} {status_text}</div>
            </div>
            <div class="content">
                <table class="meta-table">
                    <tr>
                        <td class="label">Date & Time</td>
                        <td class="value">{timestamp}</td>
                    </tr>
                    <tr>
                        <td class="label">Output Directory</td>
                        <td class="value"><code>{outdir}</code></td>
                    </tr>
                    <tr>
                        <td class="label">Samples Processed</td>
                        <td class="value">{len(samples)} sample(s)</td>
                    </tr>
                </table>

                {links_html}
                {error_html}

                <div style="margin-top: 20px;">
                    <h4 style="margin: 0 0 8px 0; color: #374151; font-size: 13px;">Target Samples:</h4>
                    <ul style="margin: 0; padding-left: 20px; font-size: 13px; color: #4b5563;">
                        {sample_list_html}
                    </ul>
                </div>
            </div>
            <div class="footer">
                DNAseq Analysis Toolkit for Cardiovascular Disease Research • Automated Snakemake Notification
            </div>
        </div>
    </body>
    </html>
    """
    return html

def send_email_notification(config, status="SUCCESS", log_file=None, samples=None, error_details=None):
    """
    Send an email notification based on workflow execution status.
    Reads configuration from the 'email' section in config.yml or CLI override 'email_recipient'.
    """
    email_cfg = config.get("email", {})
    if isinstance(email_cfg, dict):
        email_cfg = dict(email_cfg)
    else:
        email_cfg = {}

    # Check CLI override
    cli_recipient = config.get("email_recipient")
    if cli_recipient:
        email_cfg["enabled"] = True
        email_cfg["to"] = cli_recipient

    if not email_cfg or not email_cfg.get("enabled", False):
        return

    # Check notification filters
    if status.upper() == "SUCCESS" and not email_cfg.get("notify_on_success", True):
        return
    if status.upper() != "SUCCESS" and not email_cfg.get("notify_on_error", True):
        return

    to_addr = email_cfg.get("to")
    if not to_addr:
        logger.warning("Email notification enabled but no recipient ('to') address specified.")
        return

    if isinstance(to_addr, str):
        recipients = [a.strip() for a in to_addr.split(",") if a.strip()]
    else:
        recipients = list(to_addr)

    from_addr = email_cfg.get("from", "wes-pipeline@localhost")
    smtp_host = email_cfg.get("smtp_host", "localhost")
    smtp_port = int(email_cfg.get("smtp_port", 25))
    use_tls = email_cfg.get("use_tls", False)
    use_ssl = email_cfg.get("use_ssl", False)
    username = email_cfg.get("username", "")
    password = email_cfg.get("password", "")
    pipeline_name = config.get("pipeline_name", "WES Analysis Pipeline")
    outdir = config.get("outdir", "")
    samples = samples or []

    subject = f"[{status.upper()}] {pipeline_name} - {os.path.basename(outdir) or 'Run'}"
    html_content = build_html_report(status, pipeline_name, outdir, samples, log_file, error_details)
    plain_content = f"{pipeline_name} Status: {status}\nDate: {datetime.now()}\nOutput: {outdir}\nSamples: {', '.join(samples)}\nLog: {log_file or 'N/A'}\n"

    msg = MIMEMultipart("alternative")
    msg["Subject"] = subject
    msg["From"] = from_addr
    msg["To"] = ", ".join(recipients)

    msg.attach(MIMEText(plain_content, "plain"))
    msg.attach(MIMEText(html_content, "html"))

    try:
        if use_ssl:
            server = smtplib.SMTP_SSL(smtp_host, smtp_port, timeout=15)
        else:
            server = smtplib.SMTP(smtp_host, smtp_port, timeout=15)
            if use_tls:
                server.starttls()

        if username and password:
            server.login(username, password)

        server.sendmail(from_addr, recipients, msg.as_string())
        server.quit()
        print(f"\n[INFO] Email notification successfully sent to {', '.join(recipients)}")
    except Exception as e:
        print(f"\n[WARNING] Could not send email notification: {e}")

if __name__ == "__main__":
    test_config = {
        "email": {
            "enabled": True,
            "to": "test@example.com",
            "from": "wes@localhost",
            "smtp_host": "localhost",
            "smtp_port": 25,
            "notify_on_success": True,
            "notify_on_error": True,
        },
        "pipeline_name": "WES Cardio Pipeline",
        "outdir": "/mnt/bucket/Snakemake_Analysis/test_run"
    }
    send_email_notification(test_config, status="SUCCESS", samples=["134095", "137842"])
