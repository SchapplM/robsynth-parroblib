% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3RPRR1G1P1A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x8]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:23
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RPRR1G1P1A0_coriolisvec_para_pf_regmin(xP, xDP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G1P1A0_coriolisvec_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRR1G1P1A0_coriolisvec_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G1P1A0_coriolisvec_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G1P1A0_coriolisvec_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G1P1A0_coriolisvec_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G1P1A0_coriolisvec_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:23:19
% EndTime: 2020-03-09 21:23:20
% DurationCPUTime: 0.90s
% Computational Cost: add. (10688->159), mult. (7566->271), div. (1086->10), fcn. (4824->32), ass. (0->148)
t482 = pkin(7) + qJ(3,3);
t490 = sin(qJ(3,3));
t444 = pkin(1) * sin(t482) + t490 * pkin(2);
t435 = 0.1e1 / t444;
t483 = pkin(7) + qJ(3,2);
t491 = sin(qJ(3,2));
t445 = pkin(1) * sin(t483) + t491 * pkin(2);
t438 = 0.1e1 / t445;
t484 = pkin(7) + qJ(3,1);
t492 = sin(qJ(3,1));
t446 = pkin(1) * sin(t484) + t492 * pkin(2);
t441 = 0.1e1 / t446;
t436 = 0.1e1 / t444 ^ 2;
t439 = 0.1e1 / t445 ^ 2;
t442 = 0.1e1 / t446 ^ 2;
t485 = sin(pkin(7));
t563 = pkin(1) * t485;
t486 = cos(pkin(7));
t562 = pkin(2) * t486;
t493 = cos(qJ(3,3));
t561 = t493 * pkin(2);
t494 = cos(qJ(3,2));
t560 = t494 * pkin(2);
t495 = cos(qJ(3,1));
t559 = t495 * pkin(2);
t478 = legFrame(3,3) + qJ(1,3);
t465 = pkin(7) + t478;
t462 = qJ(3,3) + t465;
t450 = sin(t462);
t408 = -pkin(1) * sin(t478) - pkin(2) * sin(t465) - pkin(3) * t450;
t453 = cos(t462);
t411 = -pkin(1) * cos(t478) - pkin(2) * cos(t465) - pkin(3) * t453;
t496 = xDP(2);
t497 = xDP(1);
t403 = t408 * t496 + t411 * t497;
t499 = 0.1e1 / pkin(3);
t400 = t403 * t499 * t435;
t417 = t450 * t496 + t453 * t497;
t405 = t417 * t435;
t394 = t400 + t405;
t475 = cos(t482);
t384 = -pkin(3) * t394 + (-pkin(1) * t475 - t561) * t405;
t551 = t403 * t436;
t388 = t394 * t551;
t474 = t486 * pkin(1) + pkin(2);
t534 = t485 * t490;
t429 = -pkin(1) * t534 + t493 * t474;
t432 = t490 * t474 + t493 * t563;
t516 = t394 * (pkin(3) + t429) / t432 * t400;
t546 = t417 * t436;
t390 = t405 + t400 / 0.2e1;
t500 = pkin(1) ^ 2;
t481 = pkin(2) ^ 2 + t500;
t498 = pkin(3) ^ 2;
t530 = 0.2e1 * pkin(1);
t531 = 0.2e1 * pkin(2) * pkin(3);
t554 = (t390 * t493 * t531 + t481 * t405 + t394 * t498 + (pkin(3) * t390 * t475 + t405 * t562) * t530) * t499;
t375 = -t516 + t388 + (-t384 - t554) * t546;
t558 = t375 * t435;
t480 = legFrame(1,3) + qJ(1,1);
t467 = pkin(7) + t480;
t464 = qJ(3,1) + t467;
t452 = sin(t464);
t410 = -pkin(1) * sin(t480) - pkin(2) * sin(t467) - pkin(3) * t452;
t455 = cos(t464);
t413 = -pkin(1) * cos(t480) - pkin(2) * cos(t467) - pkin(3) * t455;
t404 = t410 * t496 + t413 * t497;
t401 = t404 * t499 * t441;
t419 = t452 * t496 + t455 * t497;
t407 = t419 * t441;
t398 = t407 + t401;
t477 = cos(t484);
t385 = -t398 * pkin(3) + (-pkin(1) * t477 - t559) * t407;
t550 = t404 * t442;
t389 = t398 * t550;
t532 = t485 * t492;
t431 = -pkin(1) * t532 + t495 * t474;
t434 = t492 * t474 + t495 * t563;
t515 = t398 * (pkin(3) + t431) / t434 * t401;
t544 = t419 * t442;
t392 = t407 + t401 / 0.2e1;
t555 = (t392 * t495 * t531 + t481 * t407 + t398 * t498 + (pkin(3) * t392 * t477 + t407 * t562) * t530) * t499;
t376 = -t515 + t389 + (-t385 - t555) * t544;
t557 = t376 * t441;
t479 = legFrame(2,3) + qJ(1,2);
t466 = pkin(7) + t479;
t463 = qJ(3,2) + t466;
t451 = sin(t463);
t409 = -pkin(1) * sin(t479) - pkin(2) * sin(t466) - pkin(3) * t451;
t454 = cos(t463);
t412 = -pkin(1) * cos(t479) - pkin(2) * cos(t466) - pkin(3) * t454;
t402 = t409 * t496 + t412 * t497;
t399 = t402 * t499 * t438;
t418 = t451 * t496 + t454 * t497;
t406 = t418 * t438;
t393 = t399 + t406;
t476 = cos(t483);
t386 = -pkin(3) * t393 + (-pkin(1) * t476 - t560) * t406;
t552 = t402 * t439;
t387 = t393 * t552;
t533 = t485 * t491;
t430 = -pkin(1) * t533 + t494 * t474;
t433 = t491 * t474 + t494 * t563;
t517 = t393 * (pkin(3) + t430) / t433 * t399;
t545 = t418 * t439;
t391 = t406 + t399 / 0.2e1;
t553 = (t391 * t494 * t531 + t481 * t406 + t393 * t498 + (pkin(3) * t391 * t476 + t406 * t562) * t530) * t499;
t377 = -t517 + t387 + (-t386 - t553) * t545;
t556 = t377 * t438;
t549 = t417 ^ 2 * t435 * t436;
t548 = t418 ^ 2 * t438 * t439;
t547 = t419 ^ 2 * t441 * t442;
t543 = t432 * t435;
t542 = t433 * t438;
t541 = t434 * t441;
t540 = t435 * t450;
t539 = t435 * t453;
t538 = t438 * t451;
t537 = t438 * t454;
t536 = t441 * t452;
t535 = t441 * t455;
t420 = -t561 + (-t486 * t493 + t534) * pkin(1);
t529 = t435 * (-t516 + 0.2e1 * t388 + (-0.2e1 * t384 - t554) * t546) * t420;
t421 = -t560 + (-t486 * t494 + t533) * pkin(1);
t528 = t438 * (-t517 + 0.2e1 * t387 + (-0.2e1 * t386 - t553) * t545) * t421;
t422 = -t559 + (-t486 * t495 + t532) * pkin(1);
t527 = t441 * (-t515 + 0.2e1 * t389 + (-0.2e1 * t385 - t555) * t544) * t422;
t526 = -0.2e1 * (t388 - t516 / 0.2e1 + (-t384 - t554 / 0.2e1) * t546) * t543;
t525 = -0.2e1 * (t387 - t517 / 0.2e1 + (-t386 - t553 / 0.2e1) * t545) * t542;
t524 = -0.2e1 * (t389 - t515 / 0.2e1 + (-t385 - t555 / 0.2e1) * t544) * t541;
t523 = (0.2e1 * t405 + t400) * t420 * t551;
t522 = (0.2e1 * t406 + t399) * t421 * t552;
t521 = (0.2e1 * t407 + t401) * t422 * t550;
t520 = -0.2e1 * t390 * t432 * t551;
t519 = -0.2e1 * t391 * t433 * t552;
t518 = -0.2e1 * t392 * t434 * t550;
t381 = -t384 * t546 + t388;
t514 = -t381 * t543 + t429 * t549;
t513 = t435 * t429 * t381 + t432 * t549;
t382 = -t385 * t544 + t389;
t512 = -t382 * t541 + t431 * t547;
t511 = t441 * t431 * t382 + t434 * t547;
t383 = -t386 * t545 + t387;
t510 = -t383 * t542 + t430 * t548;
t509 = t438 * t430 * t383 + t433 * t548;
t508 = t381 * t540 + t382 * t536 + t383 * t538;
t507 = t381 * t539 + t382 * t535 + t383 * t537;
t1 = [t507, 0, 0, t507 * t500, t375 * t539 + t376 * t535 + t377 * t537 + (t411 * t558 + t412 * t556 + t413 * t557) * t499, -t455 * t527 - t453 * t529 - t454 * t528 + (t513 * t411 + t509 * t412 + t511 * t413 + t453 * t520 + t454 * t519 + t455 * t518) * t499, t453 * t526 + t455 * t524 + t454 * t525 + (t514 * t411 + t510 * t412 + t512 * t413 + t453 * t523 + t454 * t522 + t455 * t521) * t499, 0; t508, 0, 0, t508 * t500, t375 * t540 + t376 * t536 + t377 * t538 + (t408 * t558 + t409 * t556 + t410 * t557) * t499, -t452 * t527 - t450 * t529 - t451 * t528 + (t513 * t408 + t509 * t409 + t511 * t410 + t450 * t520 + t451 * t519 + t452 * t518) * t499, t450 * t526 + t452 * t524 + t451 * t525 + (t514 * t408 + t510 * t409 + t512 * t410 + t450 * t523 + t451 * t522 + t452 * t521) * t499, 0; 0, 0, 0, 0, 0, 0, 0, 0;];
tau_reg  = t1;
