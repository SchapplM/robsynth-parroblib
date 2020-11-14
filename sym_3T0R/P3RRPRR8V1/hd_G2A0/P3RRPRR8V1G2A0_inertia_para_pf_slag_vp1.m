% Calculate inertia matrix for parallel robot
% P3RRPRR8V1G2A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d4,theta3]';
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
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 19:59
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRPRR8V1G2A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(5,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V1G2A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V1G2A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRPRR8V1G2A0_inertia_para_pf_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR8V1G2A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR8V1G2A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRPRR8V1G2A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V1G2A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V1G2A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:59:04
% EndTime: 2020-08-06 19:59:05
% DurationCPUTime: 1.08s
% Computational Cost: add. (3243->235), mult. (5142->374), div. (294->9), fcn. (2688->44), ass. (0->192)
t571 = m(2) / 0.2e1;
t570 = Icges(2,2) / 0.2e1;
t569 = Icges(3,2) / 0.2e1;
t503 = pkin(1) ^ 2;
t568 = t503 / 0.2e1;
t463 = qJ(2,3) + pkin(5);
t447 = cos(t463);
t486 = cos(qJ(2,3));
t456 = t486 * pkin(1);
t567 = pkin(2) * t447 + t456;
t464 = qJ(2,2) + pkin(5);
t448 = cos(t464);
t488 = cos(qJ(2,2));
t457 = t488 * pkin(1);
t566 = pkin(2) * t448 + t457;
t465 = qJ(2,1) + pkin(5);
t449 = cos(t465);
t490 = cos(qJ(2,1));
t458 = t490 * pkin(1);
t565 = pkin(2) * t449 + t458;
t564 = m(2) * rSges(2,3);
t563 = m(3) * rSges(3,1);
t500 = rSges(2,2) ^ 2;
t502 = rSges(2,1) ^ 2;
t562 = m(3) * t568 + (-t500 + t502) * t571 + t570 - Icges(2,1) / 0.2e1;
t499 = rSges(3,2) ^ 2;
t501 = rSges(3,1) ^ 2;
t561 = (-t499 + t501) * m(3) / 0.2e1 + t569 - Icges(3,1) / 0.2e1;
t418 = t456 + t447 * rSges(3,1) - sin(t463) * rSges(3,2);
t560 = m(3) * t418;
t419 = t457 + t448 * rSges(3,1) - sin(t464) * rSges(3,2);
t559 = m(3) * t419;
t420 = t458 + t449 * rSges(3,1) - sin(t465) * rSges(3,2);
t558 = m(3) * t420;
t474 = pkin(4) + qJ(3,3);
t466 = 0.1e1 / t474;
t557 = m(3) * t466;
t475 = pkin(4) + qJ(3,2);
t467 = 0.1e1 / t475;
t556 = m(3) * t467;
t476 = pkin(4) + qJ(3,1);
t468 = 0.1e1 / t476;
t555 = m(3) * t468;
t471 = qJ(3,3) + rSges(3,3);
t554 = m(3) * t471;
t472 = qJ(3,2) + rSges(3,3);
t553 = m(3) * t472;
t473 = qJ(3,1) + rSges(3,3);
t552 = m(3) * t473;
t469 = sin(pkin(5));
t548 = pkin(2) * t469;
t480 = sin(qJ(2,3));
t481 = sin(qJ(1,3));
t477 = legFrame(3,2);
t450 = sin(t477);
t521 = t450 * t548;
t453 = cos(t477);
t524 = t453 * t548;
t470 = cos(pkin(5));
t441 = t470 * pkin(2) + pkin(1);
t529 = t450 * t441;
t532 = t441 * t453;
t402 = (-t481 * t529 + t524) * t486 + t480 * (t481 * t521 + t532);
t507 = t441 * t486 - t480 * t548;
t421 = 0.1e1 / t507;
t547 = t402 * t421;
t482 = sin(qJ(2,2));
t483 = sin(qJ(1,2));
t478 = legFrame(2,2);
t451 = sin(t478);
t520 = t451 * t548;
t454 = cos(t478);
t523 = t454 * t548;
t528 = t451 * t441;
t531 = t441 * t454;
t403 = (-t483 * t528 + t523) * t488 + t482 * (t483 * t520 + t531);
t506 = t441 * t488 - t482 * t548;
t422 = 0.1e1 / t506;
t546 = t403 * t422;
t484 = sin(qJ(2,1));
t485 = sin(qJ(1,1));
t479 = legFrame(1,2);
t452 = sin(t479);
t519 = t452 * t548;
t455 = cos(t479);
t522 = t455 * t548;
t527 = t452 * t441;
t530 = t441 * t455;
t404 = (-t485 * t527 + t522) * t490 + t484 * (t485 * t519 + t530);
t505 = t441 * t490 - t484 * t548;
t423 = 0.1e1 / t505;
t545 = t404 * t423;
t405 = (t481 * t532 + t521) * t486 + (-t481 * t524 + t529) * t480;
t544 = t405 * t421;
t406 = (t483 * t531 + t520) * t488 + (-t483 * t523 + t528) * t482;
t543 = t406 * t422;
t407 = (t485 * t530 + t519) * t490 + (-t485 * t522 + t527) * t484;
t542 = t407 * t423;
t541 = t418 * t421;
t540 = t419 * t422;
t539 = t420 * t423;
t427 = 0.1e1 / t567;
t538 = t427 * t450;
t537 = t427 * t453;
t428 = 0.1e1 / t566;
t536 = t428 * t451;
t535 = t428 * t454;
t429 = 0.1e1 / t565;
t534 = t429 * t452;
t533 = t429 * t455;
t526 = rSges(2,2) * t564 - Icges(2,6);
t525 = t500 + t502;
t434 = rSges(3,2) * t554 - Icges(3,6);
t437 = rSges(3,1) * t554 - Icges(3,5);
t509 = -rSges(2,1) * t564 + Icges(2,5);
t393 = (-pkin(1) * t554 + t434 * t469 - t437 * t470 + t509) * t480 - t486 * (t434 * t470 + t437 * t469 + t526);
t518 = t393 * t421 * t466;
t435 = rSges(3,2) * t553 - Icges(3,6);
t438 = rSges(3,1) * t553 - Icges(3,5);
t394 = (-pkin(1) * t553 + t435 * t469 - t438 * t470 + t509) * t482 - t488 * (t435 * t470 + t438 * t469 + t526);
t517 = t394 * t422 * t467;
t436 = rSges(3,2) * t552 - Icges(3,6);
t439 = rSges(3,1) * t552 - Icges(3,5);
t395 = (-pkin(1) * t552 + t436 * t469 - t439 * t470 + t509) * t484 - t490 * (t436 * t470 + t439 * t469 + t526);
t516 = t395 * t423 * t468;
t515 = t393 * t538;
t514 = t394 * t536;
t513 = t395 * t534;
t512 = t393 * t537;
t511 = t394 * t535;
t510 = t395 * t533;
t508 = t499 / 0.2e1 + t501 / 0.2e1 + t568;
t440 = t470 * pkin(1) * t563;
t504 = Icges(1,3) + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + (0.2e1 * rSges(2,3) ^ 2 + t525) * t571 + t440 + t569 + t570 + Icges(3,1) / 0.2e1 + Icges(2,1) / 0.2e1;
t498 = 0.2e1 * qJ(2,1);
t497 = 0.2e1 * qJ(2,2);
t496 = 0.2e1 * qJ(2,3);
t491 = cos(qJ(1,1));
t489 = cos(qJ(1,2));
t487 = cos(qJ(1,3));
t462 = pkin(5) + t498;
t461 = pkin(5) + t497;
t460 = pkin(5) + t496;
t446 = -m(2) * rSges(2,1) * rSges(2,2) + Icges(2,4);
t445 = -rSges(3,2) * t563 + Icges(3,4);
t444 = 0.2e1 * t465;
t443 = 0.2e1 * t464;
t442 = 0.2e1 * t463;
t426 = t484 * t441 + t490 * t548;
t425 = t482 * t441 + t488 * t548;
t424 = t480 * t441 + t486 * t548;
t417 = t485 * t476 + t565 * t491;
t416 = t483 * t475 + t566 * t489;
t415 = t481 * t474 + t567 * t487;
t414 = -t491 * t476 + t505 * t485;
t413 = -t489 * t475 + t506 * t483;
t412 = -t487 * t474 + t507 * t481;
t411 = 0.2e1 * t440 + t525 * m(2) + Icges(2,3) + Icges(3,3) + (-0.2e1 * t469 * rSges(3,2) * pkin(1) + t499 + t501 + t503) * m(3);
t410 = (-t420 * t491 + t417) * t555;
t409 = (-t419 * t489 + t416) * t556;
t408 = (-t418 * t487 + t415) * t557;
t401 = -t414 * t452 + t426 * t455;
t400 = t414 * t455 + t426 * t452;
t399 = -t413 * t451 + t425 * t454;
t398 = t413 * t454 + t425 * t451;
t397 = -t412 * t450 + t424 * t453;
t396 = t412 * t453 + t424 * t450;
t392 = cos(t444) * t561 + cos(t498) * t562 + t445 * sin(t444) + t446 * sin(t498) + (t473 ^ 2 + (rSges(3,1) * cos(t462) + (-sin(t462) - t469) * rSges(3,2)) * pkin(1) + t508) * m(3) + t504;
t391 = cos(t443) * t561 + cos(t497) * t562 + t445 * sin(t443) + t446 * sin(t497) + (t472 ^ 2 + (rSges(3,1) * cos(t461) + (-sin(t461) - t469) * rSges(3,2)) * pkin(1) + t508) * m(3) + t504;
t390 = cos(t442) * t561 + cos(t496) * t562 + t445 * sin(t442) + t446 * sin(t496) + (t471 ^ 2 + (rSges(3,1) * cos(t460) + (-sin(t460) - t469) * rSges(3,2)) * pkin(1) + t508) * m(3) + t504;
t389 = (-t407 * t539 + t400) * t555;
t388 = (-t406 * t540 + t398) * t556;
t387 = (-t405 * t541 + t396) * t557;
t386 = (-t404 * t539 + t401) * t555;
t385 = (-t403 * t540 + t399) * t556;
t384 = (-t402 * t541 + t397) * t557;
t383 = (t392 * t491 - t417 * t558) * t468;
t382 = (t391 * t489 - t416 * t559) * t467;
t381 = (t390 * t487 - t415 * t560) * t466;
t380 = t407 * t516 + t411 * t534;
t379 = t406 * t517 + t411 * t536;
t378 = t405 * t518 + t411 * t538;
t377 = t404 * t516 + t411 * t533;
t376 = t403 * t517 + t411 * t535;
t375 = t402 * t518 + t411 * t537;
t374 = t513 + (t392 * t542 - t400 * t558) * t468;
t373 = t514 + (t391 * t543 - t398 * t559) * t467;
t372 = t515 + (t390 * t544 - t396 * t560) * t466;
t371 = t510 + (t392 * t545 - t401 * t558) * t468;
t370 = t511 + (t391 * t546 - t399 * t559) * t467;
t369 = t512 + (t390 * t547 - t397 * t560) * t466;
t1 = [t378 * t538 + t379 * t536 + t380 * t534 + m(4) + (t374 * t542 + t389 * t400) * t468 + (t373 * t543 + t388 * t398) * t467 + (t372 * t544 + t387 * t396) * t466, t378 * t537 + t379 * t535 + t380 * t533 + (t374 * t545 + t389 * t401) * t468 + (t373 * t546 + t388 * t399) * t467 + (t372 * t547 + t387 * t397) * t466, (t374 * t491 + t389 * t417) * t468 + (t373 * t489 + t388 * t416) * t467 + (t372 * t487 + t387 * t415) * t466; t375 * t538 + t376 * t536 + t377 * t534 + (t371 * t542 + t386 * t400) * t468 + (t370 * t543 + t385 * t398) * t467 + (t369 * t544 + t384 * t396) * t466, t375 * t537 + t376 * t535 + t377 * t533 + m(4) + (t371 * t545 + t386 * t401) * t468 + (t370 * t546 + t385 * t399) * t467 + (t369 * t547 + t384 * t397) * t466, (t371 * t491 + t386 * t417) * t468 + (t370 * t489 + t385 * t416) * t467 + (t369 * t487 + t384 * t415) * t466; (t383 * t542 + t400 * t410 + t491 * t513) * t468 + (t382 * t543 + t398 * t409 + t489 * t514) * t467 + (t381 * t544 + t396 * t408 + t487 * t515) * t466, (t383 * t545 + t401 * t410 + t491 * t510) * t468 + (t382 * t546 + t399 * t409 + t489 * t511) * t467 + (t381 * t547 + t397 * t408 + t487 * t512) * t466, m(4) + (t383 * t491 + t410 * t417) * t468 + (t382 * t489 + t409 * t416) * t467 + (t381 * t487 + t408 * t415) * t466;];
MX  = t1;
