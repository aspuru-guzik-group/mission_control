import abc

class BaseJobSubmissionBuilder(abc.ABCMeta):
    """Abstract Base Class for SubmissionBuilders."""

    def __init__(self, job=None, cfg=None, output_dir=None, **kwargs):
        """
        Args:
            job [dict]: job dict.
            cfg [dict]: cfg dict.
            output_dir [str]: directory in which to put submission files.
        """
        self.job = job or {}
        self.cfg = cfg or {}
        self.output_dir = output_dir

    @classmethod
    def build_job_submission(cls, job=None, cfg=None, output_dir=None,
                             **kwargs):
        """
        Args:
            job (dict): job dict.
            cfg (dict, optional): cfg dict.
            output_dir (str, optional): directory in which to put submission
                files.

        Returns:
            submission_meta (dict): dict of submission metadata.
        """
        builder = cls(job=job, cfg=cfg, output_dir=output_dir, **kwargs)
        return builder._build_job_submission()

    @abc.abstractmethod
    def _build_job_submission(self): raise NotImplementedError
    
