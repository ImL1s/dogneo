"""LLM routing layer — selects optimal backend based on task and availability.

Priority: CLI (free) → Local (offline) → Cloud API (highest quality).
Supports automatic fallback when higher-priority backends are unavailable.
"""

from __future__ import annotations

import logging
from enum import Enum

from dogneo.config import LLMConfig
from dogneo.llm.backends import (
    AnthropicBackend,
    ClaudeCodeCLIBackend,
    CodexCLIBackend,
    GeminiAPIBackend,
    GeminiCLIBackend,
    LLMBackend,
    LocalLlamaBackend,
    OpenAIBackend,
)

logger = logging.getLogger(__name__)


class TaskType(Enum):
    """Types of tasks for LLM routing."""

    SUMMARIZE = "summarize"     # Summarize analysis results
    ANALYZE = "analyze"         # Deep analysis of candidates
    TRANSLATE = "translate"     # Translate report text
    REPORT = "report"           # Generate full report
    SIMPLE = "simple"           # Simple text generation


class LLMRouter:
    """Routes LLM requests to the optimal backend.

    Three-tier priority system:
    1. CLI backends (default) — free, no API key needed
    2. Local models — offline, privacy-sensitive workloads
    3. Cloud API — complex tasks requiring highest quality

    The router automatically falls back to the next tier if
    the preferred tier is unavailable.

    Args:
        config: LLM configuration with model selections.
        force_tier: Override tier selection (for testing).
    """

    def __init__(
        self,
        config: LLMConfig | None = None,
        force_tier: str | None = None,
    ):
        self.config = config or LLMConfig()
        self.force_tier = force_tier
        self._backends: dict[str, list[LLMBackend]] = {
            "cli": [],
            "local": [],
            "cloud": [],
        }
        self._initialized = False

    def _init_backends(self) -> None:
        """Initialize all configured backends."""
        if self._initialized:
            return

        # CLI backends
        self._backends["cli"] = [
            GeminiCLIBackend(model=self.config.gemini_cli_model, timeout=self.config.cli_timeout),
            ClaudeCodeCLIBackend(model=self.config.claude_cli_model, timeout=self.config.cli_timeout),
            CodexCLIBackend(model=self.config.codex_cli_model, timeout=self.config.cli_timeout),
        ]

        # Local backend (if model path configured)
        if self.config.local_model_path:
            self._backends["local"] = [
                LocalLlamaBackend(
                    model_path=self.config.local_model_path,
                    n_ctx=self.config.local_n_ctx,
                    n_gpu_layers=self.config.local_n_gpu_layers,
                ),
            ]

        # Cloud backends (if API keys available)
        cloud_backends: list[LLMBackend] = []
        if self.config.openai_api_key:
            cloud_backends.append(OpenAIBackend(
                model=self.config.cloud_model,
                api_key=self.config.openai_api_key,
            ))
        if self.config.anthropic_api_key:
            cloud_backends.append(AnthropicBackend(
                api_key=self.config.anthropic_api_key,
            ))
        if self.config.google_api_key:
            cloud_backends.append(GeminiAPIBackend(
                api_key=self.config.google_api_key,
            ))
        self._backends["cloud"] = cloud_backends

        self._initialized = True

    def _get_tier_order(self, task_type: TaskType) -> list[str]:
        """Determine tier priority based on task type.

        For most tasks: CLI → Local → Cloud
        For complex analysis/reports: Cloud → CLI → Local
        """
        if self.force_tier:
            return [self.force_tier]

        default_tier = self.config.default_tier

        if task_type in (TaskType.ANALYZE, TaskType.REPORT):
            # Complex tasks benefit from cloud models
            if default_tier == "cloud":
                return ["cloud", "cli", "local"]

        # Default: follow configured priority
        if default_tier == "cli":
            return ["cli", "local", "cloud"]
        elif default_tier == "local":
            return ["local", "cli", "cloud"]
        elif default_tier == "cloud":
            return ["cloud", "cli", "local"]

        return ["cli", "local", "cloud"]

    def generate(
        self,
        prompt: str,
        task_type: TaskType = TaskType.SIMPLE,
        **kwargs,
    ) -> str:
        """Generate text using the best available backend.

        Automatically tries backends in tier priority order,
        falling back on failure.

        Args:
            prompt: Input prompt text.
            task_type: Type of task for routing.
            **kwargs: Backend-specific options.

        Returns:
            Generated text.

        Raises:
            RuntimeError: If all backends fail.
        """
        self._init_backends()
        tier_order = self._get_tier_order(task_type)

        errors: list[str] = []

        for tier in tier_order:
            backends = self._backends.get(tier, [])
            for backend in backends:
                try:
                    logger.debug("Trying %s backend: %s", tier, backend.name)
                    result = backend.generate(prompt, **kwargs)
                    logger.info("Generated with %s", backend.name)
                    return result
                except Exception as e:
                    error_msg = f"{backend.name}: {e}"
                    logger.debug("Backend failed: %s", error_msg)
                    errors.append(error_msg)
                    continue

        error_summary = "; ".join(errors) if errors else "No backends available"
        raise RuntimeError(f"All LLM backends failed: {error_summary}")

    def get_available_backends(self) -> dict[str, list[str]]:
        """List available backends by tier.

        Returns:
            Dict mapping tier name to list of backend names.
        """
        self._init_backends()
        result: dict[str, list[str]] = {}
        for tier, backends in self._backends.items():
            available = [b.name for b in backends if b.is_available()]
            if available:
                result[tier] = available
        return result


def create_router(config: LLMConfig | None = None) -> LLMRouter:
    """Create and return an LLM router with the given configuration.

    Args:
        config: LLM configuration. Uses defaults if None.

    Returns:
        Configured LLMRouter instance.
    """
    return LLMRouter(config=config)
