% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3RPR1G1P1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% qJ [2x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x10]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:58
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RPR1G1P1A0_coriolisvec_para_pf_regmin(xP, xDP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(2,3),zeros(3,3),zeros(3,3),zeros(4,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPR1G1P1A0_coriolisvec_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPR1G1P1A0_coriolisvec_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3RPR1G1P1A0_coriolisvec_para_pf_regmin: qJ has to be [2x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RPR1G1P1A0_coriolisvec_para_pf_regmin: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPR1G1P1A0_coriolisvec_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPR1G1P1A0_coriolisvec_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:58:12
% EndTime: 2019-05-03 14:58:13
% DurationCPUTime: 1.41s
% Computational Cost: add. (11097->207), mult. (17549->371), div. (567->9), fcn. (8218->14), ass. (0->153)
t488 = pkin(1) + pkin(2);
t502 = koppelP(3,2);
t505 = koppelP(3,1);
t458 = -qJ(2,3) * t502 + t488 * t505;
t461 = qJ(2,3) * t505 + t488 * t502;
t486 = xDP(2);
t464 = t486 * t488;
t489 = xP(3);
t472 = sin(t489);
t473 = cos(t489);
t485 = xDP(3);
t487 = xDP(1);
t410 = qJ(2,3) * t487 + t464 + (t458 * t473 - t461 * t472) * t485;
t465 = t487 * t488;
t413 = qJ(2,3) * t486 - t465 + (t458 * t472 + t461 * t473) * t485;
t476 = legFrame(3,3);
t466 = sin(t476);
t469 = cos(t476);
t479 = sin(qJ(1,3));
t482 = cos(qJ(1,3));
t395 = (t410 * t479 - t413 * t482) * t469 + (t410 * t482 + t413 * t479) * t466;
t440 = t472 * t505 + t473 * t502;
t428 = t440 * t485 - t487;
t443 = -t472 * t502 + t473 * t505;
t431 = t443 * t485 + t486;
t404 = (-t428 * t482 + t479 * t431) * t469 + t466 * (t479 * t428 + t431 * t482);
t452 = -t482 * qJ(2,3) + t479 * t488;
t455 = t479 * qJ(2,3) + t488 * t482;
t416 = t452 * t469 + t466 * t455;
t419 = -t466 * t452 + t455 * t469;
t474 = t485 ^ 2;
t490 = qJ(2,3) ^ 2;
t491 = 0.1e1 / qJ(2,3);
t492 = 0.1e1 / qJ(2,3) ^ 2;
t509 = (pkin(1) ^ 2);
t531 = -t509 + (-2 * pkin(1) - pkin(2)) * pkin(2);
t493 = t491 / t490;
t537 = t493 * t395;
t515 = -(((-t490 + t531) * t404 + t395 * t488) * t491 * t492 + t488 * t537) * t404 - (t416 * t440 + t419 * t443) * t474 * t491;
t503 = koppelP(2,2);
t506 = koppelP(2,1);
t459 = -qJ(2,2) * t503 + t488 * t506;
t462 = qJ(2,2) * t506 + t488 * t503;
t411 = qJ(2,2) * t487 + t464 + (t459 * t473 - t462 * t472) * t485;
t414 = qJ(2,2) * t486 - t465 + (t459 * t472 + t462 * t473) * t485;
t477 = legFrame(2,3);
t467 = sin(t477);
t470 = cos(t477);
t480 = sin(qJ(1,2));
t483 = cos(qJ(1,2));
t396 = (t411 * t480 - t414 * t483) * t470 + (t411 * t483 + t414 * t480) * t467;
t441 = t472 * t506 + t473 * t503;
t429 = t441 * t485 - t487;
t444 = -t472 * t503 + t473 * t506;
t432 = t444 * t485 + t486;
t405 = (-t429 * t483 + t480 * t432) * t470 + t467 * (t480 * t429 + t432 * t483);
t453 = -t483 * qJ(2,2) + t480 * t488;
t456 = t480 * qJ(2,2) + t488 * t483;
t417 = t453 * t470 + t467 * t456;
t420 = -t467 * t453 + t456 * t470;
t494 = qJ(2,2) ^ 2;
t495 = 0.1e1 / qJ(2,2);
t496 = 0.1e1 / qJ(2,2) ^ 2;
t497 = t495 / t494;
t536 = t497 * t396;
t514 = -(((-t494 + t531) * t405 + t396 * t488) * t495 * t496 + t488 * t536) * t405 - (t417 * t441 + t420 * t444) * t474 * t495;
t504 = koppelP(1,2);
t507 = koppelP(1,1);
t460 = -qJ(2,1) * t504 + t488 * t507;
t463 = qJ(2,1) * t507 + t488 * t504;
t412 = qJ(2,1) * t487 + t464 + (t460 * t473 - t463 * t472) * t485;
t415 = qJ(2,1) * t486 - t465 + (t460 * t472 + t463 * t473) * t485;
t478 = legFrame(1,3);
t468 = sin(t478);
t471 = cos(t478);
t481 = sin(qJ(1,1));
t484 = cos(qJ(1,1));
t397 = (t412 * t481 - t415 * t484) * t471 + (t412 * t484 + t415 * t481) * t468;
t442 = t472 * t507 + t473 * t504;
t430 = t442 * t485 - t487;
t445 = -t472 * t504 + t473 * t507;
t433 = t445 * t485 + t486;
t406 = (-t430 * t484 + t481 * t433) * t471 + t468 * (t481 * t430 + t433 * t484);
t454 = -t484 * qJ(2,1) + t481 * t488;
t457 = t481 * qJ(2,1) + t488 * t484;
t418 = t454 * t471 + t468 * t457;
t421 = -t468 * t454 + t457 * t471;
t498 = qJ(2,1) ^ 2;
t499 = 0.1e1 / qJ(2,1);
t500 = 0.1e1 / qJ(2,1) ^ 2;
t501 = t499 / t498;
t535 = t501 * t397;
t513 = -(((-t498 + t531) * t406 + t397 * t488) * t499 * t500 + t488 * t535) * t406 - (t418 * t442 + t421 * t445) * t474 * t499;
t434 = t466 * t482 + t469 * t479;
t435 = -t466 * t479 + t469 * t482;
t534 = t404 * t537;
t545 = t404 * t492;
t389 = -t534 + (-(-t404 * t488 + t395) * t545 + (-t434 * t440 - t435 * t443) * t474) * t491;
t554 = pkin(1) * t389;
t436 = t467 * t483 + t470 * t480;
t437 = -t467 * t480 + t470 * t483;
t533 = t405 * t536;
t544 = t405 * t496;
t390 = -t533 + (-(-t405 * t488 + t396) * t544 + (-t436 * t441 - t437 * t444) * t474) * t495;
t553 = pkin(1) * t390;
t438 = t468 * t484 + t471 * t481;
t439 = -t468 * t481 + t471 * t484;
t532 = t406 * t535;
t543 = t406 * t500;
t391 = -t532 + (-(-t406 * t488 + t397) * t543 + (-t438 * t442 - t439 * t445) * t474) * t499;
t552 = pkin(1) * t391;
t551 = t389 * t491;
t550 = t390 * t495;
t549 = t391 * t499;
t401 = t404 ^ 2;
t548 = t401 * t493;
t402 = t405 ^ 2;
t547 = t402 * t497;
t403 = t406 ^ 2;
t546 = t403 * t501;
t542 = t472 * t474;
t541 = t473 * t474;
t530 = (t515 - t554) * t491 - t401 * t492;
t529 = (t514 - t553) * t495 - t402 * t496;
t528 = (t513 - t552) * t499 - t403 * t500;
t524 = ((t490 + t509) * t389 - pkin(1) * t515) * t491 + 0.2e1 * t395 * t545;
t523 = ((t494 + t509) * t390 - pkin(1) * t514) * t495 + 0.2e1 * t396 * t544;
t522 = ((t498 + t509) * t391 - pkin(1) * t513) * t499 + 0.2e1 * t397 * t543;
t518 = 0.2e1 * t389 + 0.2e1 * t534;
t517 = 0.2e1 * t390 + 0.2e1 * t533;
t516 = 0.2e1 * t391 + 0.2e1 * t532;
t451 = t481 * t504 + t484 * t507;
t450 = t481 * t507 - t484 * t504;
t449 = t480 * t503 + t483 * t506;
t448 = t480 * t506 - t483 * t503;
t447 = t479 * t502 + t482 * t505;
t446 = t479 * t505 - t482 * t502;
t427 = t481 * t460 - t463 * t484;
t426 = t480 * t459 - t462 * t483;
t425 = t479 * t458 - t461 * t482;
t424 = t460 * t484 + t463 * t481;
t423 = t459 * t483 + t462 * t480;
t422 = t458 * t482 + t461 * t479;
t409 = (t450 * t473 - t472 * t451) * t471 + t468 * (t472 * t450 + t451 * t473);
t408 = (t448 * t473 - t472 * t449) * t470 + t467 * (t472 * t448 + t449 * t473);
t407 = (t446 * t473 - t472 * t447) * t469 + t466 * (t472 * t446 + t447 * t473);
t400 = (-t472 * t424 + t427 * t473) * t471 + (t424 * t473 + t427 * t472) * t468;
t399 = (-t472 * t423 + t426 * t473) * t470 + (t423 * t473 + t426 * t472) * t467;
t398 = (-t472 * t422 + t425 * t473) * t469 + (t422 * t473 + t425 * t472) * t466;
t387 = -t513 + 0.2e1 * t552;
t385 = -t514 + 0.2e1 * t553;
t383 = -t515 + 0.2e1 * t554;
t1 = [t435 * t551 + t437 * t550 + t439 * t549, 0, 0, (t387 * t439 - t391 * t421) * t499 + (t385 * t437 - t390 * t420) * t495 + (t383 * t435 - t389 * t419) * t491, -t419 * t548 - t420 * t547 - t421 * t546 + t435 * t518 + t437 * t517 + t439 * t516, t530 * t419 + t529 * t420 + t528 * t421 + t524 * t435 + t523 * t437 + t522 * t439, 0, -t541, t542, 0; t434 * t551 + t436 * t550 + t438 * t549, 0, 0, (t387 * t438 - t391 * t418) * t499 + (t385 * t436 - t390 * t417) * t495 + (t383 * t434 - t389 * t416) * t491, -t416 * t548 - t417 * t547 - t418 * t546 + t434 * t518 + t436 * t517 + t438 * t516, t530 * t416 + t529 * t417 + t528 * t418 + t524 * t434 + t523 * t436 + t522 * t438, 0, -t542, -t541, 0; t407 * t551 + t408 * t550 + t409 * t549, 0, 0, (t387 * t409 - t391 * t400) * t499 + (t385 * t408 - t390 * t399) * t495 + (t383 * t407 - t389 * t398) * t491, -t398 * t548 - t399 * t547 - t400 * t546 + t407 * t518 + t408 * t517 + t409 * t516, t530 * t398 + t529 * t399 + t528 * t400 + t524 * t407 + t523 * t408 + t522 * t409, 0, 0, 0, 0;];
tau_reg  = t1;
