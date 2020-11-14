% Calculate vector of centrifugal and coriolis load on the joints for
% P3PRRRR1G3A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:02
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PRRRR1G3A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G3A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR1G3A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G3A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G3A0_coriolisvec_para_pf_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR1G3A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRRR1G3A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRRR1G3A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G3A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G3A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:02:22
% EndTime: 2020-03-09 21:02:22
% DurationCPUTime: 0.64s
% Computational Cost: add. (771->130), mult. (2034->297), div. (1914->14), fcn. (2427->18), ass. (0->113)
t472 = -2 * Ifges(3,4);
t414 = cos(qJ(3,3));
t391 = 0.1e1 / t414;
t416 = cos(qJ(3,2));
t395 = 0.1e1 / t416;
t418 = cos(qJ(3,1));
t399 = 0.1e1 / t418;
t390 = t414 ^ 2;
t471 = 0.2e1 * t390;
t394 = t416 ^ 2;
t470 = 0.2e1 * t394;
t398 = t418 ^ 2;
t469 = 0.2e1 * t398;
t404 = mrSges(2,2) - mrSges(3,3);
t468 = t404 / 0.2e1;
t467 = -Ifges(3,1) - Ifges(2,3);
t408 = sin(qJ(3,3));
t466 = mrSges(3,1) * t408;
t410 = sin(qJ(3,2));
t465 = mrSges(3,1) * t410;
t412 = sin(qJ(3,1));
t464 = mrSges(3,1) * t412;
t405 = legFrame(3,2);
t380 = sin(t405);
t383 = cos(t405);
t409 = sin(qJ(2,3));
t387 = 0.1e1 / t409;
t420 = xDP(3);
t421 = xDP(2);
t422 = xDP(1);
t423 = 0.1e1 / pkin(2);
t392 = 0.1e1 / t414 ^ 2;
t415 = cos(qJ(2,3));
t443 = t392 * t408 * t415;
t361 = (-t420 * t443 + (t380 * t421 - t383 * t422) * t391) * t423 * t387;
t463 = t361 * t415;
t406 = legFrame(2,2);
t381 = sin(t406);
t384 = cos(t406);
t411 = sin(qJ(2,2));
t388 = 0.1e1 / t411;
t396 = 0.1e1 / t416 ^ 2;
t417 = cos(qJ(2,2));
t441 = t396 * t410 * t417;
t362 = (-t420 * t441 + (t381 * t421 - t384 * t422) * t395) * t423 * t388;
t462 = t362 * t417;
t407 = legFrame(1,2);
t382 = sin(t407);
t385 = cos(t407);
t413 = sin(qJ(2,1));
t389 = 0.1e1 / t413;
t400 = 0.1e1 / t418 ^ 2;
t419 = cos(qJ(2,1));
t439 = t400 * t412 * t419;
t363 = (-t420 * t439 + (t382 * t421 - t385 * t422) * t399) * t423 * t389;
t461 = t363 * t419;
t460 = (mrSges(3,2) * t414 + t466) * t409;
t459 = (mrSges(3,2) * t416 + t465) * t411;
t458 = (mrSges(3,2) * t418 + t464) * t413;
t402 = t420 ^ 2;
t457 = t402 * t423;
t424 = 0.1e1 / pkin(2) ^ 2;
t456 = t402 * t424;
t455 = t408 * t414;
t454 = t410 * t416;
t453 = t412 * t418;
t452 = t420 * t423;
t451 = t420 * t424;
t393 = t391 / t390;
t444 = t391 * t452;
t346 = ((t390 * pkin(2) * t463 - t408 * t409 * t420) * t423 * t392 * t361 + (-t409 * t361 * t455 + t415 * t444) * t393 * t452) * t387;
t358 = t361 ^ 2;
t352 = (-pkin(2) * t358 * t414 - t393 * t457) * t387;
t433 = -t414 * mrSges(3,1) + t408 * mrSges(3,2);
t364 = -(mrSges(2,1) - t433) * t415 + t409 * t404;
t379 = mrSges(3,2) * t452;
t386 = -m(1) - m(2) - m(3);
t438 = t393 * t408 * t456;
t450 = t364 * t346 + t386 * t352 - t438 * t460 + (-mrSges(2,1) * t358 + t433 * (t392 * t456 + t358)) * t409 - 0.2e1 * (t361 * t468 + t444 * t466 + t379) * t463;
t397 = t395 / t394;
t442 = t395 * t452;
t347 = ((t394 * pkin(2) * t462 - t410 * t411 * t420) * t423 * t396 * t362 + (-t411 * t362 * t454 + t417 * t442) * t397 * t452) * t388;
t359 = t362 ^ 2;
t353 = (-pkin(2) * t359 * t416 - t397 * t457) * t388;
t432 = -t416 * mrSges(3,1) + t410 * mrSges(3,2);
t365 = -(mrSges(2,1) - t432) * t417 + t411 * t404;
t437 = t397 * t410 * t456;
t449 = t365 * t347 + t386 * t353 - t437 * t459 + (-mrSges(2,1) * t359 + t432 * (t396 * t456 + t359)) * t411 - 0.2e1 * (t362 * t468 + t442 * t465 + t379) * t462;
t401 = t399 / t398;
t440 = t399 * t452;
t348 = ((t398 * pkin(2) * t461 - t412 * t413 * t420) * t423 * t400 * t363 + (-t413 * t363 * t453 + t419 * t440) * t401 * t452) * t389;
t360 = t363 ^ 2;
t354 = (-pkin(2) * t360 * t418 - t401 * t457) * t389;
t431 = -t418 * mrSges(3,1) + t412 * mrSges(3,2);
t366 = -(mrSges(2,1) - t431) * t419 + t413 * t404;
t436 = t401 * t412 * t456;
t448 = t366 * t348 + t386 * t354 - t436 * t458 + (-mrSges(2,1) * t360 + t431 * (t400 * t456 + t360)) * t413 - 0.2e1 * (t363 * t468 + t440 * t464 + t379) * t461;
t403 = Ifges(3,1) - Ifges(3,2);
t434 = -Ifges(3,6) * t452 / 0.2e1;
t435 = Ifges(3,5) * t452 / 0.2e1;
t447 = ((t361 * t403 * t408 + t391 * t435) * t414 + t391 * t408 * t434 + (t471 - 0.1e1) * t361 * Ifges(3,4)) * t451;
t446 = ((t362 * t403 * t410 + t395 * t435) * t416 + t395 * t410 * t434 + (t470 - 0.1e1) * t362 * Ifges(3,4)) * t451;
t445 = ((t363 * t403 * t412 + t399 * t435) * t418 + t399 * t412 * t434 + (t469 - 0.1e1) * t363 * Ifges(3,4)) * t451;
t373 = -Ifges(3,5) * t408 - Ifges(3,6) * t414;
t337 = t364 * t352 + (t403 * t390 + t455 * t472 + t467) * t346 - t373 * t438;
t430 = t337 * t391 * t423 + 0.2e1 * t392 * t447;
t374 = -Ifges(3,5) * t410 - Ifges(3,6) * t416;
t338 = t365 * t353 + (t403 * t394 + t454 * t472 + t467) * t347 - t374 * t437;
t429 = t338 * t395 * t423 + 0.2e1 * t396 * t446;
t375 = -Ifges(3,5) * t412 - Ifges(3,6) * t418;
t339 = t366 * t354 + (t403 * t398 + t453 * t472 + t467) * t348 - t375 * t436;
t428 = t339 * t399 * t423 + 0.2e1 * t400 * t445;
t1 = [(-t428 * t385 + t448 * (t413 * t382 + t385 * t419)) * t389 + (-t429 * t384 + t449 * (t411 * t381 + t384 * t417)) * t388 + (-t430 * t383 + t450 * (t409 * t380 + t383 * t415)) * t387; (t428 * t382 + t448 * (-t382 * t419 + t385 * t413)) * t389 + (t429 * t381 + t449 * (-t381 * t417 + t384 * t411)) * t388 + (t430 * t380 + t450 * (-t380 * t415 + t383 * t409)) * t387; (-0.2e1 * t401 * t419 * t445 + t448 * t399) * t412 * t389 + (-0.2e1 * t397 * t417 * t446 + t449 * t395) * t410 * t388 + (-0.2e1 * t393 * t415 * t447 + t450 * t391) * t408 * t387 + (-t387 * t337 * t443 - t388 * t338 * t441 - t389 * t339 * t439 + (Ifges(3,3) * t436 + t354 * t458 + t375 * t348 - (Ifges(3,4) * t469 + t403 * t453 - Ifges(3,4)) * t360) * t399 + (Ifges(3,3) * t437 + t353 * t459 + t374 * t347 - (Ifges(3,4) * t470 + t403 * t454 - Ifges(3,4)) * t359) * t395 + (Ifges(3,3) * t438 + t352 * t460 + t373 * t346 - (Ifges(3,4) * t471 + t403 * t455 - Ifges(3,4)) * t358) * t391) * t423;];
taucX  = t1;
