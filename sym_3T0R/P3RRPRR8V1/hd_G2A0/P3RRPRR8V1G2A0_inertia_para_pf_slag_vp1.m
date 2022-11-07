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
% Datum: 2022-11-04 17:05
% Revision: e482436b586c4f286726c907c195760c5ac72455 (2022-11-03)
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
% StartTime: 2022-11-04 17:04:14
% EndTime: 2022-11-04 17:04:15
% DurationCPUTime: 1.01s
% Computational Cost: add. (3243->235), mult. (5142->374), div. (294->9), fcn. (2688->44), ass. (0->192)
t562 = m(2) / 0.2e1;
t561 = Icges(2,2) / 0.2e1;
t560 = Icges(3,2) / 0.2e1;
t494 = pkin(1) ^ 2;
t559 = t494 / 0.2e1;
t453 = qJ(2,3) + pkin(5);
t438 = cos(t453);
t477 = cos(qJ(2,3));
t447 = t477 * pkin(1);
t558 = pkin(2) * t438 + t447;
t455 = qJ(2,2) + pkin(5);
t439 = cos(t455);
t479 = cos(qJ(2,2));
t448 = t479 * pkin(1);
t557 = pkin(2) * t439 + t448;
t456 = qJ(2,1) + pkin(5);
t440 = cos(t456);
t481 = cos(qJ(2,1));
t449 = t481 * pkin(1);
t556 = pkin(2) * t440 + t449;
t555 = m(2) * rSges(2,3);
t554 = m(3) * rSges(3,1);
t491 = rSges(2,2) ^ 2;
t493 = rSges(2,1) ^ 2;
t553 = m(3) * t559 + (-t491 + t493) * t562 + t561 - Icges(2,1) / 0.2e1;
t490 = rSges(3,2) ^ 2;
t492 = rSges(3,1) ^ 2;
t552 = (-t490 + t492) * m(3) / 0.2e1 + t560 - Icges(3,1) / 0.2e1;
t409 = t447 + t438 * rSges(3,1) - sin(t453) * rSges(3,2);
t551 = m(3) * t409;
t410 = t448 + t439 * rSges(3,1) - sin(t455) * rSges(3,2);
t550 = m(3) * t410;
t411 = t449 + t440 * rSges(3,1) - sin(t456) * rSges(3,2);
t549 = m(3) * t411;
t465 = pkin(4) + qJ(3,3);
t457 = 0.1e1 / t465;
t548 = m(3) * t457;
t466 = pkin(4) + qJ(3,2);
t458 = 0.1e1 / t466;
t547 = m(3) * t458;
t467 = pkin(4) + qJ(3,1);
t459 = 0.1e1 / t467;
t546 = m(3) * t459;
t462 = rSges(3,3) + qJ(3,3);
t545 = m(3) * t462;
t463 = rSges(3,3) + qJ(3,2);
t544 = m(3) * t463;
t464 = rSges(3,3) + qJ(3,1);
t543 = m(3) * t464;
t460 = sin(pkin(5));
t539 = pkin(2) * t460;
t471 = sin(qJ(2,3));
t472 = sin(qJ(1,3));
t468 = legFrame(3,2);
t441 = sin(t468);
t512 = t441 * t539;
t444 = cos(t468);
t515 = t444 * t539;
t461 = cos(pkin(5));
t432 = t461 * pkin(2) + pkin(1);
t520 = t441 * t432;
t523 = t432 * t444;
t393 = (-t472 * t520 + t515) * t477 + (t472 * t512 + t523) * t471;
t498 = t432 * t477 - t471 * t539;
t412 = 0.1e1 / t498;
t538 = t393 * t412;
t473 = sin(qJ(2,2));
t474 = sin(qJ(1,2));
t469 = legFrame(2,2);
t442 = sin(t469);
t511 = t442 * t539;
t445 = cos(t469);
t514 = t445 * t539;
t519 = t442 * t432;
t522 = t432 * t445;
t394 = (-t474 * t519 + t514) * t479 + (t474 * t511 + t522) * t473;
t497 = t432 * t479 - t473 * t539;
t413 = 0.1e1 / t497;
t537 = t394 * t413;
t475 = sin(qJ(2,1));
t476 = sin(qJ(1,1));
t470 = legFrame(1,2);
t443 = sin(t470);
t510 = t443 * t539;
t446 = cos(t470);
t513 = t446 * t539;
t518 = t443 * t432;
t521 = t432 * t446;
t395 = (-t476 * t518 + t513) * t481 + (t476 * t510 + t521) * t475;
t496 = t432 * t481 - t475 * t539;
t414 = 0.1e1 / t496;
t536 = t395 * t414;
t396 = (t472 * t523 + t512) * t477 + (-t472 * t515 + t520) * t471;
t535 = t396 * t412;
t397 = (t474 * t522 + t511) * t479 + (-t474 * t514 + t519) * t473;
t534 = t397 * t413;
t398 = (t476 * t521 + t510) * t481 + (-t476 * t513 + t518) * t475;
t533 = t398 * t414;
t532 = t409 * t412;
t531 = t410 * t413;
t530 = t411 * t414;
t418 = 0.1e1 / t558;
t529 = t418 * t441;
t528 = t418 * t444;
t419 = 0.1e1 / t557;
t527 = t419 * t442;
t526 = t419 * t445;
t420 = 0.1e1 / t556;
t525 = t420 * t443;
t524 = t420 * t446;
t517 = rSges(2,2) * t555 - Icges(2,6);
t516 = t491 + t493;
t425 = rSges(3,2) * t545 - Icges(3,6);
t428 = rSges(3,1) * t545 - Icges(3,5);
t500 = -rSges(2,1) * t555 + Icges(2,5);
t384 = (-pkin(1) * t545 + t425 * t460 - t428 * t461 + t500) * t471 - (t425 * t461 + t428 * t460 + t517) * t477;
t509 = t384 * t412 * t457;
t426 = rSges(3,2) * t544 - Icges(3,6);
t429 = rSges(3,1) * t544 - Icges(3,5);
t385 = (-pkin(1) * t544 + t426 * t460 - t429 * t461 + t500) * t473 - (t426 * t461 + t429 * t460 + t517) * t479;
t508 = t385 * t413 * t458;
t427 = rSges(3,2) * t543 - Icges(3,6);
t430 = rSges(3,1) * t543 - Icges(3,5);
t386 = (-pkin(1) * t543 + t427 * t460 - t430 * t461 + t500) * t475 - (t427 * t461 + t430 * t460 + t517) * t481;
t507 = t386 * t414 * t459;
t506 = t384 * t529;
t505 = t385 * t527;
t504 = t386 * t525;
t503 = t384 * t528;
t502 = t385 * t526;
t501 = t386 * t524;
t499 = t490 / 0.2e1 + t492 / 0.2e1 + t559;
t431 = pkin(1) * t461 * t554;
t495 = Icges(1,3) + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + (0.2e1 * rSges(2,3) ^ 2 + t516) * t562 + t431 + t560 + t561 + Icges(3,1) / 0.2e1 + Icges(2,1) / 0.2e1;
t489 = 0.2e1 * qJ(2,1);
t488 = 0.2e1 * qJ(2,2);
t487 = 0.2e1 * qJ(2,3);
t482 = cos(qJ(1,1));
t480 = cos(qJ(1,2));
t478 = cos(qJ(1,3));
t454 = t488 + pkin(5);
t452 = pkin(5) + t489;
t451 = pkin(5) + t487;
t437 = -m(2) * rSges(2,1) * rSges(2,2) + Icges(2,4);
t436 = -rSges(3,2) * t554 + Icges(3,4);
t435 = 0.2e1 * t456;
t434 = 0.2e1 * t455;
t433 = 0.2e1 * t453;
t417 = t475 * t432 + t481 * t539;
t416 = t473 * t432 + t479 * t539;
t415 = t471 * t432 + t477 * t539;
t408 = t476 * t467 + t556 * t482;
t407 = t474 * t466 + t557 * t480;
t406 = t472 * t465 + t558 * t478;
t405 = -t482 * t467 + t496 * t476;
t404 = -t480 * t466 + t497 * t474;
t403 = -t478 * t465 + t498 * t472;
t402 = 0.2e1 * t431 + t516 * m(2) + Icges(2,3) + Icges(3,3) + (-0.2e1 * pkin(1) * rSges(3,2) * t460 + t490 + t492 + t494) * m(3);
t401 = (-t411 * t482 + t408) * t546;
t400 = (-t410 * t480 + t407) * t547;
t399 = (-t409 * t478 + t406) * t548;
t392 = -t405 * t443 + t446 * t417;
t391 = t405 * t446 + t443 * t417;
t390 = -t404 * t442 + t445 * t416;
t389 = t404 * t445 + t442 * t416;
t388 = -t403 * t441 + t444 * t415;
t387 = t403 * t444 + t441 * t415;
t383 = cos(t435) * t552 + cos(t489) * t553 + t436 * sin(t435) + t437 * sin(t489) + (t464 ^ 2 + (rSges(3,1) * cos(t452) + (-sin(t452) - t460) * rSges(3,2)) * pkin(1) + t499) * m(3) + t495;
t382 = cos(t434) * t552 + cos(t488) * t553 + t436 * sin(t434) + t437 * sin(t488) + (t463 ^ 2 + (rSges(3,1) * cos(t454) + (-sin(t454) - t460) * rSges(3,2)) * pkin(1) + t499) * m(3) + t495;
t381 = cos(t433) * t552 + cos(t487) * t553 + t436 * sin(t433) + t437 * sin(t487) + (t462 ^ 2 + (rSges(3,1) * cos(t451) + (-sin(t451) - t460) * rSges(3,2)) * pkin(1) + t499) * m(3) + t495;
t380 = (-t398 * t530 + t391) * t546;
t379 = (-t397 * t531 + t389) * t547;
t378 = (-t396 * t532 + t387) * t548;
t377 = (-t395 * t530 + t392) * t546;
t376 = (-t394 * t531 + t390) * t547;
t375 = (-t393 * t532 + t388) * t548;
t374 = (t383 * t482 - t408 * t549) * t459;
t373 = (t382 * t480 - t407 * t550) * t458;
t372 = (t381 * t478 - t406 * t551) * t457;
t371 = t398 * t507 + t402 * t525;
t370 = t397 * t508 + t402 * t527;
t369 = t396 * t509 + t402 * t529;
t368 = t395 * t507 + t402 * t524;
t367 = t394 * t508 + t402 * t526;
t366 = t393 * t509 + t402 * t528;
t365 = t504 + (t383 * t533 - t391 * t549) * t459;
t364 = t505 + (t382 * t534 - t389 * t550) * t458;
t363 = t506 + (t381 * t535 - t387 * t551) * t457;
t362 = t501 + (t383 * t536 - t392 * t549) * t459;
t361 = t502 + (t382 * t537 - t390 * t550) * t458;
t360 = t503 + (t381 * t538 - t388 * t551) * t457;
t1 = [t369 * t529 + t370 * t527 + t371 * t525 + m(4) + (t365 * t533 + t380 * t391) * t459 + (t364 * t534 + t379 * t389) * t458 + (t363 * t535 + t378 * t387) * t457, t369 * t528 + t370 * t526 + t371 * t524 + (t365 * t536 + t380 * t392) * t459 + (t364 * t537 + t379 * t390) * t458 + (t363 * t538 + t378 * t388) * t457, (t365 * t482 + t380 * t408) * t459 + (t364 * t480 + t379 * t407) * t458 + (t363 * t478 + t378 * t406) * t457; t366 * t529 + t367 * t527 + t368 * t525 + (t362 * t533 + t377 * t391) * t459 + (t361 * t534 + t376 * t389) * t458 + (t360 * t535 + t375 * t387) * t457, t366 * t528 + t367 * t526 + t368 * t524 + m(4) + (t362 * t536 + t377 * t392) * t459 + (t361 * t537 + t376 * t390) * t458 + (t360 * t538 + t375 * t388) * t457, (t362 * t482 + t377 * t408) * t459 + (t361 * t480 + t376 * t407) * t458 + (t360 * t478 + t375 * t406) * t457; (t374 * t533 + t391 * t401 + t482 * t504) * t459 + (t373 * t534 + t389 * t400 + t480 * t505) * t458 + (t372 * t535 + t387 * t399 + t478 * t506) * t457, (t374 * t536 + t392 * t401 + t482 * t501) * t459 + (t373 * t537 + t390 * t400 + t480 * t502) * t458 + (t372 * t538 + t388 * t399 + t478 * t503) * t457, m(4) + (t374 * t482 + t401 * t408) * t459 + (t373 * t480 + t400 * t407) * t458 + (t372 * t478 + t399 * t406) * t457;];
MX  = t1;
