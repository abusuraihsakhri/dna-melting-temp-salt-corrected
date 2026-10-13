"""
FastAPI REST API Server for Dna Melting Temp Salt Corrected.
"""
from typing import Dict, Any, List
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
from .base import AuditLogger, PHIGuard
from .models import SystemTaskPayload, ConsensusDossier
from .supervisor import SystemSupervisor
from mg2_correction import calculate_tm_with_mg2

supervisor = SystemSupervisor(model_provider="mock")

app = FastAPI(
    title="Dna Melting Temp Salt Corrected API",
    description="Enterprise Distributed Component Platform (Clinical & Biomedical AI)",
    version="3.0.0-ENTERPRISE",
)


class ChatRequest(BaseModel):
    query: str


@app.get("/health")
def health():
    return {"status": "HEALTHY", "service": "dna-melting-temp-salt-corrected", "domain": "Clinical & Biomedical AI", "standard": "CAP / CLSI / ISO Standards", "version": "3.0.0-ENTERPRISE"}


@app.get("/metrics")
def metrics():
    return {
        "dossiers_processed_total": len(supervisor.dossier_registry),
        "audit_blocks_total": len(AuditLogger.get_trail()),
        "system_status": "NOMINAL_OPTIMAL"
    }


class TmRequest(BaseModel):
    sequence: str
    na_mM: float = 50.0
    mg_mM: float = 1.5
    oligo_conc_nM: float = 250.0


@app.post("/api/tm")
def api_tm(req: TmRequest):
    """Estimate salt-corrected Tm for a perfectly matched DNA duplex."""
    result = calculate_tm_with_mg2(
        req.sequence, na_mM=req.na_mM, mg_mM=req.mg_mM,
        oligo_conc_nM=req.oligo_conc_nM,
    )
    if "error" in result:
        raise HTTPException(status_code=422, detail=result["error"])
    return result


@app.post("/api/audit")
def api_audit(payload: SystemTaskPayload):
    dossier = supervisor.process_task(payload)
    return dossier.to_dict()


@app.post("/api/chat")
def api_chat(req: ChatRequest):
    try:
        ans = supervisor.query_supervisory_chat(req.query)
        return {"response": ans}
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@app.get("/api/audit/logs")
def api_audit_logs():
    return {"audit_trail": AuditLogger.get_trail(), "verified": AuditLogger.verify_integrity()}
