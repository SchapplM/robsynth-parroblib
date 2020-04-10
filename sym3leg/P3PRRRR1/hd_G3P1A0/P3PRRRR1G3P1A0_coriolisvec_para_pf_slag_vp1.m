% Calculate vector of centrifugal and coriolis load on the joints for
% P3PRRRR1G3P1A0
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
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
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

function taucX = P3PRRRR1G3P1A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G3P1A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR1G3P1A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G3P1A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G3P1A0_coriolisvec_para_pf_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR1G3P1A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRRR1G3P1A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRRR1G3P1A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G3P1A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G3P1A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:02:17
% EndTime: 2020-03-09 21:02:18
% DurationCPUTime: 0.81s
% Computational Cost: add. (948->142), mult. (2586->313), div. (1914->14), fcn. (2418->24), ass. (0->124)
t450 = sin(qJ(3,1));
t456 = cos(qJ(3,1));
t522 = (rSges(3,1) * t456 - rSges(3,2) * t450) * m(3);
t448 = sin(qJ(3,2));
t454 = cos(qJ(3,2));
t521 = (rSges(3,1) * t454 - rSges(3,2) * t448) * m(3);
t446 = sin(qJ(3,3));
t452 = cos(qJ(3,3));
t520 = (rSges(3,1) * t452 - rSges(3,2) * t446) * m(3);
t431 = 0.1e1 / t452;
t435 = 0.1e1 / t454;
t439 = 0.1e1 / t456;
t419 = m(3) * rSges(3,1) * rSges(3,2) - Icges(3,4);
t519 = 0.2e1 * t419;
t460 = m(2) * rSges(2,1);
t518 = rSges(3,3) * m(3);
t467 = rSges(3,2) ^ 2;
t468 = rSges(3,1) ^ 2;
t412 = (-t467 + t468) * m(3) - Icges(3,1) + Icges(3,2);
t517 = -t412 / 0.2e1;
t417 = rSges(3,2) * t518 - Icges(3,6);
t516 = -t417 / 0.4e1;
t418 = rSges(3,1) * t518 - Icges(3,5);
t515 = t418 / 0.4e1;
t443 = legFrame(3,2);
t420 = sin(t443);
t423 = cos(t443);
t447 = sin(qJ(2,3));
t427 = 0.1e1 / t447;
t461 = xDP(3);
t462 = xDP(2);
t463 = xDP(1);
t469 = 0.1e1 / pkin(2);
t432 = 0.1e1 / t452 ^ 2;
t453 = cos(qJ(2,3));
t496 = t432 * t446 * t453;
t393 = (-t461 * t496 + (t420 * t462 - t423 * t463) * t431) * t469 * t427;
t514 = t393 * t453;
t444 = legFrame(2,2);
t421 = sin(t444);
t424 = cos(t444);
t449 = sin(qJ(2,2));
t428 = 0.1e1 / t449;
t436 = 0.1e1 / t454 ^ 2;
t455 = cos(qJ(2,2));
t494 = t436 * t448 * t455;
t394 = (-t461 * t494 + (t421 * t462 - t424 * t463) * t435) * t469 * t428;
t513 = t394 * t455;
t445 = legFrame(1,2);
t422 = sin(t445);
t425 = cos(t445);
t451 = sin(qJ(2,1));
t429 = 0.1e1 / t451;
t440 = 0.1e1 / t456 ^ 2;
t457 = cos(qJ(2,1));
t492 = t440 * t450 * t457;
t395 = (-t461 * t492 + (t422 * t462 - t425 * t463) * t439) * t469 * t429;
t512 = t395 * t457;
t442 = t461 ^ 2;
t511 = t442 * t469;
t470 = 0.1e1 / pkin(2) ^ 2;
t510 = t442 * t470;
t509 = t446 * t452;
t508 = t448 * t454;
t507 = t450 * t456;
t506 = t461 * t469;
t430 = t452 ^ 2;
t433 = t431 / t430;
t497 = t431 * t506;
t383 = ((pkin(2) * t430 * t514 - t446 * t447 * t461) * t469 * t432 * t393 + (-t393 * t447 * t509 + t453 * t497) * t433 * t506) * t427;
t390 = t393 ^ 2;
t384 = (-pkin(2) * t390 * t452 - t433 * t511) * t427;
t416 = m(2) * rSges(2,2) - t518;
t396 = -(t460 + t520) * t453 + t447 * t416;
t488 = 0.2e1 * m(3) * t506;
t415 = rSges(3,2) * t488;
t426 = -m(1) - m(2) - m(3);
t481 = rSges(3,1) * t488;
t484 = t433 * t446 * t510;
t500 = m(3) * (rSges(3,1) * t446 + rSges(3,2) * t452) * t447;
t505 = t396 * t383 + t426 * t384 - t484 * t500 + (-t390 * t460 - (t432 * t510 + t390) * t520) * t447 - (t431 * t446 * t481 + t393 * t416 + t415) * t514;
t434 = t454 ^ 2;
t437 = t435 / t434;
t495 = t435 * t506;
t381 = ((pkin(2) * t434 * t513 - t448 * t449 * t461) * t469 * t436 * t394 + (-t394 * t449 * t508 + t455 * t495) * t437 * t506) * t428;
t391 = t394 ^ 2;
t385 = (-pkin(2) * t391 * t454 - t437 * t511) * t428;
t397 = -(t460 + t521) * t455 + t449 * t416;
t483 = t437 * t448 * t510;
t499 = m(3) * (rSges(3,1) * t448 + rSges(3,2) * t454) * t449;
t504 = t397 * t381 + t426 * t385 - t483 * t499 + (-t391 * t460 - (t436 * t510 + t391) * t521) * t449 - (t435 * t448 * t481 + t394 * t416 + t415) * t513;
t438 = t456 ^ 2;
t441 = t439 / t438;
t493 = t439 * t506;
t382 = ((pkin(2) * t438 * t512 - t450 * t451 * t461) * t469 * t440 * t395 + (-t395 * t451 * t507 + t457 * t493) * t441 * t506) * t429;
t392 = t395 ^ 2;
t386 = (-pkin(2) * t392 * t456 - t441 * t511) * t429;
t398 = -(t460 + t522) * t457 + t451 * t416;
t482 = t441 * t450 * t510;
t498 = m(3) * (rSges(3,1) * t450 + rSges(3,2) * t456) * t451;
t503 = t398 * t382 + t426 * t386 - t482 * t498 + (-t392 * t460 - (t440 * t510 + t392) * t522) * t451 - (t439 * t450 * t481 + t395 * t416 + t415) * t512;
t502 = t467 + t468;
t501 = 0.4e1 * t461 * t470;
t491 = t412 * t509;
t490 = t412 * t508;
t489 = t412 * t507;
t487 = ((t446 * t516 + t452 * t515) * t497 + (t491 / 0.2e1 + (t430 - 0.1e1 / 0.2e1) * t419) * t393) * t501;
t486 = ((t448 * t516 + t454 * t515) * t495 + (t490 / 0.2e1 + (t434 - 0.1e1 / 0.2e1) * t419) * t394) * t501;
t485 = ((t450 * t516 + t456 * t515) * t493 + (t489 / 0.2e1 + (t438 - 0.1e1 / 0.2e1) * t419) * t395) * t501;
t477 = -(rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) - (0.2e1 * rSges(3,3) ^ 2 + t502) * m(3) / 0.2e1 - Icges(3,2) / 0.2e1 - Icges(3,1) / 0.2e1 - Icges(2,3);
t399 = t417 * t452 + t418 * t446;
t464 = 0.2e1 * qJ(3,3);
t369 = t396 * t384 + (cos(t464) * t517 + t419 * sin(t464) + t477) * t383 - t399 * t484;
t476 = -t369 * t431 * t469 + t432 * t487;
t400 = t417 * t454 + t418 * t448;
t465 = 0.2e1 * qJ(3,2);
t370 = t397 * t385 + (cos(t465) * t517 + t419 * sin(t465) + t477) * t381 - t400 * t483;
t475 = -t370 * t435 * t469 + t436 * t486;
t401 = t417 * t456 + t418 * t450;
t466 = 0.2e1 * qJ(3,1);
t371 = t398 * t386 + (cos(t466) * t517 + t419 * sin(t466) + t477) * t382 - t401 * t482;
t474 = -t371 * t439 * t469 + t440 * t485;
t414 = -t502 * m(3) - Icges(3,3);
t1 = [(t474 * t425 + t503 * (t422 * t451 + t425 * t457)) * t429 + (t475 * t424 + t504 * (t421 * t449 + t424 * t455)) * t428 + (t476 * t423 + t505 * (t420 * t447 + t423 * t453)) * t427; (-t474 * t422 + t503 * (-t422 * t457 + t425 * t451)) * t429 + (-t475 * t421 + t504 * (-t421 * t455 + t424 * t449)) * t428 + (-t476 * t420 + t505 * (-t420 * t453 + t423 * t447)) * t427; (t441 * t457 * t485 + t503 * t439) * t450 * t429 + (t437 * t455 * t486 + t504 * t435) * t448 * t428 + (t433 * t453 * t487 + t505 * t431) * t446 * t427 + (-t427 * t369 * t496 - t428 * t370 * t494 - t429 * t371 * t492 + (-t414 * t482 + t386 * t498 + t401 * t382 + t392 * (t438 * t519 - t419 + t489)) * t439 + (-t414 * t483 + t385 * t499 + t400 * t381 + t391 * (t434 * t519 - t419 + t490)) * t435 + (-t414 * t484 + t384 * t500 + t399 * t383 + t390 * (t430 * t519 - t419 + t491)) * t431) * t469;];
taucX  = t1;
