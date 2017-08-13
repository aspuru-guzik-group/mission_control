import unittest

import sqlalchemy as _sqla
import sqlalchemy.orm as _sqla_orm

from mc.utils import hash_utils

from .. import models


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        self.engine = _sqla.create_engine('sqlite://')
        models.utils.Base.metadata.create_all(self.engine)
        self.Session = _sqla_orm.sessionmaker(bind=self.engine)
        self.session = self.Session()
        self.job_type = 'some.job_type'
        self.job_params = {'param_%s' % i: 'value_%s' % i for i in range(3)}
        self.props = {'str_prop': 'str', 'int_prop': 1, 'bool_prop': True}
        self.tags = {'tag_%s' % i for i in range(3)}
        self.artifact_meta = {'some': 'artifact_meta'}


class CreateJobTestCase(BaseTestCase):
    def setUp(self):
        super().setUp()
        self.job = models.Job(
            job_type=self.job_type,
            job_params=self.job_params,
            props=self.props,
            tags=self.tags,
            artifact_meta=self.artifact_meta
        )
        self.session.add(self.job)
        self.session.commit()

    def test_has_timestamps(self):
        self.assertTrue(self.job.created is not None)
        self.assertTrue(self.job.modified is not None)

    def test_key_startswith_job(self):
        self.assertTrue(self.job.key.startswith('job:'))

    def test_has_job_type(self):
        self.assertTrue(self.job.job_type, self.job_type)

    def test_has_job_params(self):
        self.assertTrue(self.job.job_params, self.job_params)

    def test_has_props(self):
        self.assertTrue(self.job.props, self.props)

    def test_has_tags(self):
        self.assertTrue(self.job.tags, self.tags)

    def test_has_default_status(self):
        self.assertTrue(self.job.status, 'PENDING')

    def test_has_artifact_meta(self):
        self.assertTrue(self.job.artifact_meta, self.artifact_meta)


class HashTestCase(BaseTestCase):
    def test_job_hash_tracks_job_type_and_job_params(self):
        job = models.Job(job_type='some_job_type',
                         job_params='some_job_params')
        expected_job_hash = hash_utils.hash_obj({
            'job_type': job.job_type,
            'job_params': job.job_params,
        })
        self.assertEqual(job.job_hash, expected_job_hash)
        self.session.add(job)
        self.session.commit()
        fetched_job = self.session.query(models.Job).first()
        self.assertEqual(fetched_job.job_hash, expected_job_hash)

        job.job_type = 'new_job_type'
        expected_job_hash = hash_utils.hash_obj({
            'job_type': job.job_type,
            'job_params': job.job_params,
        })
        self.assertEqual(job.job_hash, expected_job_hash)
        self.session.add(job)
        self.session.commit()
        fetched_job = self.session.query(models.Job).first()
        self.assertEqual(fetched_job.job_hash, expected_job_hash)

        job.job_params = 'new_job_params'
        expected_job_hash = hash_utils.hash_obj({
            'job_type': job.job_type,
            'job_params': job.job_params,
        })
        self.assertEqual(job.job_hash, expected_job_hash)
        self.session.add(job)
        self.session.commit()
        fetched_job = self.session.query(models.Job).first()
        self.assertEqual(fetched_job.job_hash, expected_job_hash)
