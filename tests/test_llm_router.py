"""Tests for dogneo.llm.router — LLM routing and tier selection."""
from __future__ import annotations

import pytest

from dogneo.config import LLMConfig
from dogneo.llm.backends import LLMBackend
from dogneo.llm.router import LLMRouter, TaskType, create_router

# ---------------------------------------------------------------------------
# Fake backend for testing
# ---------------------------------------------------------------------------

class FakeBackend(LLMBackend):
    """A controllable test backend."""

    def __init__(self, name: str, tier: str, available: bool = True, output: str = "ok"):
        self._name = name
        self._tier = tier
        self._available = available
        self._output = output

    @property
    def name(self) -> str:
        return self._name

    @property
    def tier(self) -> str:
        return self._tier

    def generate(self, prompt: str, **kwargs) -> str:
        if not self._available:
            raise RuntimeError(f"{self._name} unavailable")
        return self._output

    def is_available(self) -> bool:
        return self._available


# ---------------------------------------------------------------------------
# Router tier selection
# ---------------------------------------------------------------------------

class TestTierOrder:
    """Tests for tier priority based on task type."""

    def test_default_cli_first(self):
        router = LLMRouter(config=LLMConfig())
        order = router._get_tier_order(TaskType.SIMPLE)
        assert order[0] == "cli"

    def test_force_tier(self):
        router = LLMRouter(config=LLMConfig(), force_tier="cloud")
        order = router._get_tier_order(TaskType.SIMPLE)
        assert order == ["cloud"]

    def test_complex_task_cloud_first(self):
        cfg = LLMConfig()
        cfg.default_tier = "cloud"
        router = LLMRouter(config=cfg)
        order = router._get_tier_order(TaskType.ANALYZE)
        assert order[0] == "cloud"


# ---------------------------------------------------------------------------
# Fallback chain
# ---------------------------------------------------------------------------

class TestFallbackChain:
    """Tests for automatic backend fallback."""

    def test_uses_first_available(self):
        router = LLMRouter()
        router._initialized = True
        router._backends = {
            "cli": [FakeBackend("cli1", "cli", available=True, output="cli_out")],
            "local": [],
            "cloud": [],
        }
        result = router.generate("test prompt")
        assert result == "cli_out"

    def test_fallback_on_failure(self):
        router = LLMRouter()
        router._initialized = True
        router._backends = {
            "cli": [FakeBackend("cli1", "cli", available=False)],
            "local": [FakeBackend("local1", "local", available=True, output="local_out")],
            "cloud": [],
        }
        result = router.generate("test prompt")
        assert result == "local_out"

    def test_all_backends_fail(self):
        router = LLMRouter()
        router._initialized = True
        router._backends = {
            "cli": [FakeBackend("cli1", "cli", available=False)],
            "local": [FakeBackend("local1", "local", available=False)],
            "cloud": [FakeBackend("cloud1", "cloud", available=False)],
        }
        with pytest.raises(RuntimeError, match="All LLM backends failed"):
            router.generate("test prompt")

    def test_no_backends(self):
        router = LLMRouter()
        router._initialized = True
        router._backends = {"cli": [], "local": [], "cloud": []}
        with pytest.raises(RuntimeError, match="No backends available"):
            router.generate("test prompt")


# ---------------------------------------------------------------------------
# get_available_backends
# ---------------------------------------------------------------------------

class TestAvailableBackends:
    """Tests for backend availability listing."""

    def test_lists_available(self):
        router = LLMRouter()
        router._initialized = True
        router._backends = {
            "cli": [
                FakeBackend("gemini", "cli", available=True),
                FakeBackend("claude", "cli", available=False),
            ],
            "local": [],
            "cloud": [FakeBackend("openai", "cloud", available=True)],
        }
        avail = router.get_available_backends()
        assert "cli" in avail
        assert avail["cli"] == ["gemini"]
        assert avail["cloud"] == ["openai"]

    def test_empty_when_none_available(self):
        router = LLMRouter()
        router._initialized = True
        router._backends = {
            "cli": [FakeBackend("x", "cli", available=False)],
            "local": [],
            "cloud": [],
        }
        avail = router.get_available_backends()
        assert avail == {}


# ---------------------------------------------------------------------------
# create_router helper
# ---------------------------------------------------------------------------

class TestCreateRouter:
    def test_returns_router(self):
        router = create_router()
        assert isinstance(router, LLMRouter)

    def test_with_config(self):
        cfg = LLMConfig()
        cfg.default_tier = "local"
        router = create_router(config=cfg)
        assert router.config.default_tier == "local"
