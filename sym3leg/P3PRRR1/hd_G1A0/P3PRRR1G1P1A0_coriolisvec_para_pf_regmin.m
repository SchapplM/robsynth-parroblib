% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3PRRR1G1P1A0
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
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2020-03-09 21:15
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3PRRR1G1P1A0_coriolisvec_para_pf_regmin(xP, xDP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G1P1A0_coriolisvec_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRR1G1P1A0_coriolisvec_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G1P1A0_coriolisvec_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G1P1A0_coriolisvec_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G1P1A0_coriolisvec_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G1P1A0_coriolisvec_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:14:54
% EndTime: 2020-03-09 21:14:56
% DurationCPUTime: 1.32s
% Computational Cost: add. (14578->145), mult. (10122->295), div. (1608->9), fcn. (10128->24), ass. (0->156)
t464 = pkin(7) + qJ(2,3);
t455 = qJ(3,3) + t464;
t443 = sin(t455);
t446 = cos(t455);
t449 = sin(t464);
t452 = cos(t464);
t422 = t443 * t452 - t449 * t446;
t561 = 0.1e1 / t422;
t465 = pkin(7) + qJ(2,2);
t456 = qJ(3,2) + t465;
t444 = sin(t456);
t447 = cos(t456);
t450 = sin(t465);
t453 = cos(t465);
t423 = t444 * t453 - t450 * t447;
t560 = 0.1e1 / t423;
t466 = pkin(7) + qJ(2,1);
t457 = qJ(3,1) + t466;
t445 = sin(t457);
t448 = cos(t457);
t451 = sin(t466);
t454 = cos(t466);
t424 = t445 * t454 - t451 * t448;
t559 = 0.1e1 / t424;
t414 = 0.1e1 / t422 ^ 2;
t417 = 0.1e1 / t423 ^ 2;
t420 = 0.1e1 / t424 ^ 2;
t480 = 0.1e1 / pkin(2);
t479 = 0.1e1 / pkin(3);
t467 = legFrame(3,3);
t458 = sin(t467);
t461 = cos(t467);
t476 = xDP(2);
t477 = xDP(1);
t434 = t458 * t476 + t461 * t477;
t435 = t458 * t477 - t461 * t476;
t401 = -t434 * t446 + t435 * t443;
t549 = (-(t434 * t452 - t449 * t435) * pkin(2) + t401 * pkin(3)) * t561;
t529 = t479 * t549;
t543 = t401 * t561;
t390 = (t529 - t543) * t480;
t493 = t443 * t449 + t446 * t452;
t380 = pkin(3) * t390 - t493 * t543;
t530 = t561 * t549;
t504 = t390 * t530;
t481 = 0.1e1 / pkin(2) ^ 2;
t531 = t479 * t481;
t383 = (t493 * pkin(2) + pkin(3)) * t504 * t531;
t542 = t401 * t414;
t478 = pkin(3) ^ 2;
t506 = 0.2e1 * pkin(3) * t480;
t552 = (-t390 * t478 + (t543 - t493 * (-t543 + t529 / 0.2e1) * t506) * pkin(2)) * t479;
t371 = -t383 + (t504 + (-t380 - t552) * t542) * t481;
t558 = t371 * t561;
t468 = legFrame(2,3);
t459 = sin(t468);
t462 = cos(t468);
t436 = t459 * t476 + t462 * t477;
t437 = t459 * t477 - t462 * t476;
t402 = -t436 * t447 + t437 * t444;
t548 = (-(t436 * t453 - t450 * t437) * pkin(2) + t402 * pkin(3)) * t560;
t527 = t479 * t548;
t541 = t402 * t560;
t392 = (t527 - t541) * t480;
t491 = t444 * t450 + t447 * t453;
t381 = pkin(3) * t392 - t491 * t541;
t528 = t560 * t548;
t502 = t392 * t528;
t384 = (t491 * pkin(2) + pkin(3)) * t502 * t531;
t540 = t402 * t417;
t551 = (-t392 * t478 + (t541 - t491 * (-t541 + t527 / 0.2e1) * t506) * pkin(2)) * t479;
t372 = -t384 + (t502 + (-t381 - t551) * t540) * t481;
t557 = t372 * t560;
t469 = legFrame(1,3);
t460 = sin(t469);
t463 = cos(t469);
t438 = t460 * t476 + t463 * t477;
t439 = t460 * t477 - t463 * t476;
t403 = -t438 * t448 + t439 * t445;
t547 = (-(t438 * t454 - t451 * t439) * pkin(2) + t403 * pkin(3)) * t559;
t525 = t479 * t547;
t539 = t403 * t559;
t394 = (t525 - t539) * t480;
t489 = t445 * t451 + t448 * t454;
t382 = t394 * pkin(3) - t489 * t539;
t526 = t559 * t547;
t500 = t394 * t526;
t385 = (t489 * pkin(2) + pkin(3)) * t500 * t531;
t538 = t403 * t420;
t550 = (-t394 * t478 + (t539 - t489 * (-t539 + t525 / 0.2e1) * t506) * pkin(2)) * t479;
t373 = -t385 + (t500 + (-t382 - t550) * t538) * t481;
t556 = t373 * t559;
t374 = (-t380 * t542 + t504) * t481;
t555 = t374 * t561;
t375 = (-t381 * t540 + t502) * t481;
t554 = t375 * t560;
t376 = (-t382 * t538 + t500) * t481;
t553 = t376 * t559;
t546 = t401 ^ 2 * t561 * t414;
t545 = t402 ^ 2 * t560 * t417;
t544 = t403 ^ 2 * t559 * t420;
t492 = t461 * t443 + t458 * t446;
t537 = t561 * t492;
t429 = t443 * t458 - t461 * t446;
t536 = t561 * t429;
t490 = t462 * t444 + t459 * t447;
t535 = t560 * t490;
t431 = t444 * t459 - t462 * t447;
t534 = t560 * t431;
t488 = t463 * t445 + t460 * t448;
t533 = t559 * t488;
t433 = t445 * t460 - t463 * t448;
t532 = t559 * t433;
t470 = sin(qJ(3,3));
t524 = t470 * t546;
t473 = cos(qJ(3,3));
t523 = t473 * t546;
t471 = sin(qJ(3,2));
t522 = t471 * t545;
t474 = cos(qJ(3,2));
t521 = t474 * t545;
t472 = sin(qJ(3,1));
t520 = t472 * t544;
t475 = cos(qJ(3,1));
t519 = t475 * t544;
t518 = t470 * t555;
t517 = t473 * t555;
t516 = t471 * t554;
t515 = t474 * t554;
t514 = t472 * t553;
t513 = t475 * t553;
t368 = -t383 + (0.2e1 * t504 + (-0.2e1 * t380 - t552) * t542) * t481;
t512 = t368 * t537;
t511 = t368 * t536;
t369 = -t384 + (0.2e1 * t502 + (-0.2e1 * t381 - t551) * t540) * t481;
t510 = t369 * t535;
t509 = t369 * t534;
t370 = -t385 + (0.2e1 * t500 + (-0.2e1 * t382 - t550) * t538) * t481;
t508 = t370 * t533;
t507 = t370 * t532;
t505 = (t529 - 0.2e1 * t543) * t480 * t530;
t503 = (t527 - 0.2e1 * t541) * t480 * t528;
t501 = (t525 - 0.2e1 * t539) * t480 * t526;
t499 = t470 * t505;
t498 = t473 * t505;
t497 = t471 * t503;
t496 = t474 * t503;
t495 = t472 * t501;
t494 = t475 * t501;
t409 = -pkin(2) * (t451 * t460 - t463 * t454) - t433 * pkin(3);
t408 = -pkin(2) * (t450 * t459 - t462 * t453) - t431 * pkin(3);
t407 = -pkin(2) * (t449 * t458 - t461 * t452) - t429 * pkin(3);
t406 = pkin(2) * (t451 * t463 + t460 * t454) + t488 * pkin(3);
t405 = pkin(2) * (t450 * t462 + t459 * t453) + t490 * pkin(3);
t404 = pkin(2) * (t449 * t461 + t458 * t452) + t492 * pkin(3);
t1 = [0, (-t374 * t536 - t375 * t534 - t376 * t532) * t480, 0, 0, (-t371 * t536 - t372 * t534 - t373 * t532 + (-t407 * t558 - t408 * t557 - t409 * t556) * t479) * t480, -t473 * t511 - t474 * t509 - t475 * t507 + (-t407 * t517 - t408 * t515 - t409 * t513 + (-t407 * t524 - t408 * t522 - t409 * t520) * t481 + (t429 * t499 + t431 * t497 + t433 * t495) * t480) * t479, t470 * t511 + t471 * t509 + t472 * t507 + (t407 * t518 + t408 * t516 + t409 * t514 + (-t407 * t523 - t408 * t521 - t409 * t519) * t481 + (t429 * t498 + t431 * t496 + t433 * t494) * t480) * t479, 0; 0, (t374 * t537 + t375 * t535 + t376 * t533) * t480, 0, 0, (t371 * t537 + t372 * t535 + t373 * t533 + (-t404 * t558 - t405 * t557 - t406 * t556) * t479) * t480, t473 * t512 + t474 * t510 + t475 * t508 + (-t404 * t517 - t405 * t515 - t406 * t513 + (-t404 * t524 - t405 * t522 - t406 * t520) * t481 + (-t488 * t495 - t490 * t497 - t492 * t499) * t480) * t479, -t470 * t512 - t471 * t510 - t472 * t508 + (t404 * t518 + t405 * t516 + t406 * t514 + (-t404 * t523 - t405 * t521 - t406 * t519) * t481 + (-t488 * t494 - t490 * t496 - t492 * t498) * t480) * t479, 0; 0, 0, 0, 0, 0, 0, 0, 0;];
tau_reg  = t1;
