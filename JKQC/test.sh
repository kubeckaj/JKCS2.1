awk '
/if cfg.wandb.watch_model:/ {
    print "    best_val = {\"loss\": float(\"inf\")}"
    print "    "
    print "    @validator.on(Events.EPOCH_COMPLETED)"
    print "    def update_best_val_loss(engine):"
    print "        current = engine.state.metrics.get(\"loss\")"
    print "        if current is not None and current < best_val[\"loss\"]:"
    print "            best_val[\"loss\"] = current"
    print "            wandb.run.summary[\"best_val/loss\"] = current"
    print "    "
}
{ print }
' JKCS/AIMNET/lib/python3.11/site-packages/aimnet/train/utils.py > tmp
