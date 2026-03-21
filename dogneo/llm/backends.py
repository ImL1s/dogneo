"""LLM backend implementations.

Three tiers of backends:
1. CLI (free, no API key) — Gemini CLI, Claude Code, Codex
2. Local (offline, privacy) — llama-cpp-python GGUF
3. Cloud API (highest quality) — OpenAI, Anthropic, Gemini API
"""

from __future__ import annotations

import logging
import time
from abc import ABC, abstractmethod

from dogneo.llm.cli_wrapper import call_ai_cli, CLIResult

logger = logging.getLogger(__name__)


class LLMBackend(ABC):
    """Abstract base class for LLM backends."""

    @property
    @abstractmethod
    def name(self) -> str:
        """Backend display name."""

    @property
    @abstractmethod
    def tier(self) -> str:
        """Backend tier: 'cli', 'local', or 'cloud'."""

    @abstractmethod
    def generate(self, prompt: str, **kwargs) -> str:
        """Generate text from a prompt.

        Args:
            prompt: Input prompt text.
            **kwargs: Backend-specific options.

        Returns:
            Generated text.

        Raises:
            RuntimeError: If generation fails.
        """

    def is_available(self) -> bool:
        """Check if this backend is ready to use."""
        return True


# ---------------------------------------------------------------------------
# Tier 1: CLI Backends (free, no API key)
# ---------------------------------------------------------------------------

class GeminiCLIBackend(LLMBackend):
    """Gemini CLI backend via subprocess.

    Free tier: 60 RPM, 1000 RPD. No API key required.
    Install: npm install -g @google/gemini-cli
    """

    def __init__(self, model: str = "gemini-2.5-flash", timeout: int = 180):
        self.model = model
        self.timeout = timeout

    @property
    def name(self) -> str:
        return f"gemini-cli:{self.model}"

    @property
    def tier(self) -> str:
        return "cli"

    def generate(self, prompt: str, **kwargs) -> str:
        result = call_ai_cli(
            prompt=prompt,
            tool="gemini",
            model=self.model,
            timeout=kwargs.get("timeout", self.timeout),
            output_format=kwargs.get("output_format", "text"),
        )
        if not result.success:
            raise RuntimeError(f"Gemini CLI failed: {result.error}")
        return result.output


class ClaudeCodeCLIBackend(LLMBackend):
    """Claude Code CLI backend via subprocess.

    Requires Anthropic Max/Pro subscription.
    Install: npm install -g @anthropic-ai/claude-code
    """

    def __init__(self, model: str = "claude-sonnet-4-6", timeout: int = 180):
        self.model = model
        self.timeout = timeout

    @property
    def name(self) -> str:
        return f"claude-cli:{self.model}"

    @property
    def tier(self) -> str:
        return "cli"

    def generate(self, prompt: str, **kwargs) -> str:
        result = call_ai_cli(
            prompt=prompt,
            tool="claude",
            model=self.model,
            timeout=kwargs.get("timeout", self.timeout),
            output_format=kwargs.get("output_format", "text"),
        )
        if not result.success:
            raise RuntimeError(f"Claude CLI failed: {result.error}")
        return result.output


class CodexCLIBackend(LLMBackend):
    """Codex CLI backend via subprocess.

    Requires ChatGPT Plus/Pro subscription.
    Install: npm install -g @openai/codex
    """

    def __init__(self, model: str = "gpt-5.4", timeout: int = 180):
        self.model = model
        self.timeout = timeout

    @property
    def name(self) -> str:
        return f"codex-cli:{self.model}"

    @property
    def tier(self) -> str:
        return "cli"

    def generate(self, prompt: str, **kwargs) -> str:
        result = call_ai_cli(
            prompt=prompt,
            tool="codex",
            model=self.model,
            timeout=kwargs.get("timeout", self.timeout),
        )
        if not result.success:
            raise RuntimeError(f"Codex CLI failed: {result.error}")
        return result.output


# ---------------------------------------------------------------------------
# Tier 2: Local Backend (offline, privacy)
# ---------------------------------------------------------------------------

class LocalLlamaBackend(LLMBackend):
    """Local LLM via llama-cpp-python (GGUF format).

    Fully offline, no data leaves the machine.
    Requires: pip install llama-cpp-python
    """

    def __init__(
        self,
        model_path: str,
        n_ctx: int = 4096,
        n_gpu_layers: int = -1,
    ):
        self.model_path = model_path
        self.n_ctx = n_ctx
        self.n_gpu_layers = n_gpu_layers
        self._llm = None

    @property
    def name(self) -> str:
        return f"local:{self.model_path.split('/')[-1]}"

    @property
    def tier(self) -> str:
        return "local"

    def _ensure_loaded(self) -> None:
        """Lazy-load the model."""
        if self._llm is None:
            try:
                from llama_cpp import Llama
                self._llm = Llama(
                    model_path=self.model_path,
                    n_ctx=self.n_ctx,
                    n_gpu_layers=self.n_gpu_layers,
                    verbose=False,
                )
                logger.info("Loaded local model: %s", self.model_path)
            except ImportError:
                raise RuntimeError(
                    "llama-cpp-python not installed. Install: pip install llama-cpp-python"
                )

    def generate(self, prompt: str, **kwargs) -> str:
        self._ensure_loaded()
        max_tokens = kwargs.get("max_tokens", 2048)
        temperature = kwargs.get("temperature", 0.7)

        result = self._llm(
            prompt,
            max_tokens=max_tokens,
            temperature=temperature,
            echo=False,
        )

        return result["choices"][0]["text"].strip()

    def is_available(self) -> bool:
        try:
            from pathlib import Path
            return Path(self.model_path).exists()
        except Exception:
            return False


# ---------------------------------------------------------------------------
# Tier 3: Cloud API Backends (highest quality)
# ---------------------------------------------------------------------------

class OpenAIBackend(LLMBackend):
    """OpenAI API backend (GPT-4o, etc.).

    Requires: OPENAI_API_KEY environment variable.
    """

    def __init__(self, model: str = "gpt-4o", api_key: str = ""):
        self.model = model
        self.api_key = api_key
        self._client = None

    @property
    def name(self) -> str:
        return f"openai:{self.model}"

    @property
    def tier(self) -> str:
        return "cloud"

    def _ensure_client(self) -> None:
        if self._client is None:
            try:
                import openai
                self._client = openai.OpenAI(api_key=self.api_key or None)
            except ImportError:
                raise RuntimeError("openai package not installed. Install: pip install openai")

    def generate(self, prompt: str, **kwargs) -> str:
        self._ensure_client()
        max_retries = kwargs.get("max_retries", 3)

        for attempt in range(max_retries):
            try:
                response = self._client.chat.completions.create(
                    model=self.model,
                    messages=[{"role": "user", "content": prompt}],
                    max_tokens=kwargs.get("max_tokens", 2048),
                    temperature=kwargs.get("temperature", 0.7),
                )
                return response.choices[0].message.content.strip()
            except Exception as e:
                if "429" in str(e) and attempt < max_retries - 1:
                    wait = 2 ** (attempt + 1)
                    logger.warning("Rate limited, waiting %ds...", wait)
                    time.sleep(wait)
                    continue
                raise RuntimeError(f"OpenAI API error: {e}")

        raise RuntimeError("Max retries exceeded")

    def is_available(self) -> bool:
        import os
        return bool(self.api_key or os.environ.get("OPENAI_API_KEY"))


class AnthropicBackend(LLMBackend):
    """Anthropic Claude API backend.

    Requires: ANTHROPIC_API_KEY environment variable.
    """

    def __init__(self, model: str = "claude-sonnet-4-6", api_key: str = ""):
        self.model = model
        self.api_key = api_key
        self._client = None

    @property
    def name(self) -> str:
        return f"anthropic:{self.model}"

    @property
    def tier(self) -> str:
        return "cloud"

    def _ensure_client(self) -> None:
        if self._client is None:
            try:
                import anthropic
                self._client = anthropic.Anthropic(api_key=self.api_key or None)
            except ImportError:
                raise RuntimeError("anthropic package not installed. Install: pip install anthropic")

    def generate(self, prompt: str, **kwargs) -> str:
        self._ensure_client()
        max_retries = kwargs.get("max_retries", 3)

        for attempt in range(max_retries):
            try:
                response = self._client.messages.create(
                    model=self.model,
                    max_tokens=kwargs.get("max_tokens", 2048),
                    messages=[{"role": "user", "content": prompt}],
                )
                return response.content[0].text.strip()
            except Exception as e:
                if "429" in str(e) and attempt < max_retries - 1:
                    wait = 2 ** (attempt + 1)
                    logger.warning("Rate limited, waiting %ds...", wait)
                    time.sleep(wait)
                    continue
                raise RuntimeError(f"Anthropic API error: {e}")

        raise RuntimeError("Max retries exceeded")

    def is_available(self) -> bool:
        import os
        return bool(self.api_key or os.environ.get("ANTHROPIC_API_KEY"))


class GeminiAPIBackend(LLMBackend):
    """Google Gemini API backend.

    Requires: GOOGLE_API_KEY environment variable.
    """

    def __init__(self, model: str = "gemini-2.5-flash", api_key: str = ""):
        self.model = model
        self.api_key = api_key
        self._client = None

    @property
    def name(self) -> str:
        return f"gemini-api:{self.model}"

    @property
    def tier(self) -> str:
        return "cloud"

    def _ensure_client(self) -> None:
        if self._client is None:
            try:
                import google.generativeai as genai
                genai.configure(api_key=self.api_key or None)
                self._client = genai.GenerativeModel(self.model)
            except ImportError:
                raise RuntimeError(
                    "google-generativeai not installed. Install: pip install google-generativeai"
                )

    def generate(self, prompt: str, **kwargs) -> str:
        self._ensure_client()
        try:
            response = self._client.generate_content(prompt)
            return response.text.strip()
        except Exception as e:
            raise RuntimeError(f"Gemini API error: {e}")

    def is_available(self) -> bool:
        import os
        return bool(self.api_key or os.environ.get("GOOGLE_API_KEY"))
