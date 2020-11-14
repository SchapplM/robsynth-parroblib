% Calculate inertia matrix for parallel robot
% P3PRRRR8V1G1A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d2,d4,theta1]';
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
% Datum: 2020-08-06 16:50
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3PRRRR8V1G1A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G1A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G1A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G1A0_inertia_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V1G1A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRRR8V1G1A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRRR8V1G1A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G1A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G1A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 16:49:43
% EndTime: 2020-08-06 16:49:44
% DurationCPUTime: 1.14s
% Computational Cost: add. (2313->181), mult. (5835->376), div. (360->7), fcn. (5844->28), ass. (0->191)
t585 = m(3) / 0.2e1;
t584 = Icges(3,2) / 0.2e1;
t491 = cos(qJ(3,1));
t473 = 0.1e1 / t491;
t486 = sin(qJ(2,1));
t492 = cos(qJ(2,1));
t573 = pkin(2) * t491;
t449 = -t492 * pkin(5) + t486 * t573;
t475 = sin(pkin(3));
t477 = cos(pkin(3));
t485 = sin(qJ(3,1));
t576 = pkin(2) * t485;
t512 = 0.1e1 / (t449 * t475 + t477 * t576);
t558 = t512 * t473;
t489 = cos(qJ(3,2));
t472 = 0.1e1 / t489;
t484 = sin(qJ(2,2));
t490 = cos(qJ(2,2));
t574 = pkin(2) * t489;
t448 = -t490 * pkin(5) + t484 * t574;
t483 = sin(qJ(3,2));
t577 = pkin(2) * t483;
t513 = 0.1e1 / (t448 * t475 + t477 * t577);
t559 = t513 * t472;
t487 = cos(qJ(3,3));
t471 = 0.1e1 / t487;
t482 = sin(qJ(2,3));
t488 = cos(qJ(2,3));
t575 = pkin(2) * t487;
t447 = -t488 * pkin(5) + t482 * t575;
t481 = sin(qJ(3,3));
t578 = pkin(2) * t481;
t514 = 0.1e1 / (t447 * t475 + t477 * t578);
t560 = t514 * t471;
t583 = m(3) * rSges(3,3);
t499 = rSges(3,2) ^ 2;
t500 = rSges(3,1) ^ 2;
t582 = (-t499 + t500) * t585 - Icges(3,1) / 0.2e1 + t584;
t508 = rSges(3,1) * t487 - rSges(3,2) * t481;
t581 = m(3) * (t508 * t477 - t475 * t482 * (t481 * rSges(3,1) + t487 * rSges(3,2)));
t507 = rSges(3,1) * t489 - rSges(3,2) * t483;
t580 = m(3) * (t507 * t477 - t475 * t484 * (t483 * rSges(3,1) + t489 * rSges(3,2)));
t506 = rSges(3,1) * t491 - rSges(3,2) * t485;
t579 = m(3) * (t506 * t477 - t475 * t486 * (t485 * rSges(3,1) + t491 * rSges(3,2)));
t478 = legFrame(3,3);
t464 = sin(t478);
t467 = cos(t478);
t474 = sin(pkin(6));
t476 = cos(pkin(6));
t435 = -t474 * t464 + t467 * t476;
t548 = t477 * t482;
t441 = t474 * t488 + t476 * t548;
t444 = -t474 * t548 + t476 * t488;
t551 = t475 * t487;
t405 = (-t441 * t467 - t464 * t444) * t481 - t435 * t551;
t572 = t405 * t471;
t479 = legFrame(2,3);
t465 = sin(t479);
t468 = cos(t479);
t436 = -t474 * t465 + t468 * t476;
t547 = t477 * t484;
t442 = t474 * t490 + t476 * t547;
t445 = -t474 * t547 + t476 * t490;
t550 = t475 * t489;
t406 = (-t442 * t468 - t465 * t445) * t483 - t436 * t550;
t571 = t406 * t472;
t480 = legFrame(1,3);
t466 = sin(t480);
t469 = cos(t480);
t437 = -t474 * t466 + t469 * t476;
t546 = t477 * t486;
t443 = t474 * t492 + t476 * t546;
t446 = -t474 * t546 + t476 * t492;
t549 = t475 * t491;
t407 = (-t443 * t469 - t466 * t446) * t485 - t437 * t549;
t570 = t407 * t473;
t438 = t476 * t464 + t467 * t474;
t408 = (-t464 * t441 + t444 * t467) * t481 - t438 * t551;
t569 = t408 * t471;
t439 = t476 * t465 + t468 * t474;
t409 = (-t465 * t442 + t445 * t468) * t483 - t439 * t550;
t568 = t409 * t472;
t440 = t476 * t466 + t469 * t474;
t410 = (-t466 * t443 + t446 * t469) * t485 - t440 * t549;
t567 = t410 * t473;
t463 = -m(3) * rSges(3,1) * rSges(3,2) + Icges(3,4);
t496 = 0.2e1 * qJ(3,3);
t533 = t499 + t500;
t505 = Icges(2,3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + (0.2e1 * rSges(3,3) ^ 2 + t533) * t585 + t584 + Icges(3,1) / 0.2e1;
t566 = (cos(t496) * t582 + t463 * sin(t496) + t505) * t471;
t497 = 0.2e1 * qJ(3,2);
t565 = (cos(t497) * t582 + t463 * sin(t497) + t505) * t472;
t498 = 0.2e1 * qJ(3,1);
t564 = (cos(t498) * t582 + t463 * sin(t498) + t505) * t473;
t460 = m(2) * rSges(2,2) - t583;
t495 = m(2) * rSges(2,1);
t563 = ((t508 * m(3) + t495) * t488 - t482 * t460) * t475;
t562 = ((t507 * m(3) + t495) * t490 - t484 * t460) * t475;
t561 = ((t506 * m(3) + t495) * t492 - t486 * t460) * t475;
t470 = m(1) + m(2) + m(3);
t557 = t514 * t470;
t556 = t513 * t470;
t555 = t512 * t470;
t461 = -rSges(3,2) * t583 + Icges(3,6);
t462 = rSges(3,1) * t583 - Icges(3,5);
t432 = t461 * t487 - t462 * t481;
t554 = t432 * t471;
t433 = t461 * t489 - t462 * t483;
t553 = t433 * t472;
t434 = t461 * t491 - t462 * t485;
t552 = t434 * t473;
t545 = t482 * t435;
t544 = t482 * t438;
t543 = t484 * t436;
t542 = t484 * t439;
t541 = t486 * t437;
t540 = t486 * t440;
t539 = t488 * t435;
t538 = t488 * t438;
t537 = t490 * t436;
t536 = t490 * t439;
t535 = t492 * t437;
t534 = t492 * t440;
t399 = -(t477 * t539 - t544) * t575 - pkin(5) * (t477 * t545 + t538);
t532 = t399 * t560;
t400 = -(t477 * t538 + t545) * t575 - (t477 * t544 - t539) * pkin(5);
t531 = t400 * t560;
t401 = -(t477 * t537 - t542) * t574 - pkin(5) * (t477 * t543 + t536);
t530 = t401 * t559;
t402 = -(t477 * t536 + t543) * t574 - (t477 * t542 - t537) * pkin(5);
t529 = t402 * t559;
t403 = -(t477 * t535 - t540) * t573 - pkin(5) * (t477 * t541 + t534);
t528 = t403 * t558;
t404 = -(t477 * t534 + t541) * t573 - (t477 * t540 - t535) * pkin(5);
t527 = t404 * t558;
t501 = 0.1e1 / pkin(2);
t526 = t501 * t560;
t525 = t501 * t559;
t524 = t501 * t558;
t450 = pkin(5) * t482 + t488 * t575;
t504 = -t447 * t477 + t475 * t578;
t414 = t450 * t476 + t504 * t474;
t417 = t474 * t450 - t504 * t476;
t393 = t414 * t467 - t464 * t417;
t511 = t526 * t581;
t523 = t563 * t560;
t363 = t393 * t557 + t399 * t511 + t405 * t523;
t451 = pkin(5) * t484 + t490 * t574;
t503 = -t448 * t477 + t475 * t577;
t415 = t451 * t476 + t503 * t474;
t418 = t474 * t451 - t503 * t476;
t394 = t415 * t468 - t465 * t418;
t510 = t525 * t580;
t522 = t562 * t559;
t364 = t394 * t556 + t401 * t510 + t406 * t522;
t452 = pkin(5) * t486 + t492 * t573;
t502 = -t449 * t477 + t475 * t576;
t416 = t452 * t476 + t502 * t474;
t419 = t474 * t452 - t502 * t476;
t395 = t416 * t469 - t466 * t419;
t509 = t524 * t579;
t521 = t561 * t558;
t365 = t395 * t555 + t403 * t509 + t407 * t521;
t396 = t414 * t464 + t417 * t467;
t366 = t396 * t557 + t400 * t511 + t408 * t523;
t397 = t415 * t465 + t418 * t468;
t367 = t397 * t556 + t402 * t510 + t409 * t522;
t398 = t416 * t466 + t419 * t469;
t368 = t398 * t555 + t404 * t509 + t410 * t521;
t520 = t432 * t526;
t456 = t533 * m(3) + Icges(3,3);
t519 = t456 * t526;
t518 = t433 * t525;
t517 = t456 * t525;
t516 = t434 * t524;
t515 = t456 * t524;
t374 = t404 * t515 + (t398 * t579 + t410 * t552) * t512;
t373 = t402 * t517 + (t397 * t580 + t409 * t553) * t513;
t372 = t400 * t519 + (t396 * t581 + t408 * t554) * t514;
t371 = t403 * t515 + (t395 * t579 + t407 * t552) * t512;
t370 = t401 * t517 + (t394 * t580 + t406 * t553) * t513;
t369 = t399 * t519 + (t393 * t581 + t405 * t554) * t514;
t362 = t404 * t516 + (t398 * t561 + t410 * t564) * t512;
t361 = t402 * t518 + (t397 * t562 + t409 * t565) * t513;
t360 = t400 * t520 + (t396 * t563 + t408 * t566) * t514;
t359 = t403 * t516 + (t395 * t561 + t407 * t564) * t512;
t358 = t401 * t518 + (t394 * t562 + t406 * t565) * t513;
t357 = t399 * t520 + (t393 * t563 + t405 * t566) * t514;
t356 = t368 + t367 + t366;
t355 = t365 + t364 + t363;
t1 = [m(4) + (t359 * t570 + t365 * t395) * t512 + (t358 * t571 + t364 * t394) * t513 + (t357 * t572 + t363 * t393) * t514 + (t369 * t532 + t370 * t530 + t371 * t528) * t501, (t359 * t567 + t365 * t398) * t512 + (t358 * t568 + t364 * t397) * t513 + (t357 * t569 + t363 * t396) * t514 + (t369 * t531 + t370 * t529 + t371 * t527) * t501, t355; (t362 * t570 + t368 * t395) * t512 + (t361 * t571 + t367 * t394) * t513 + (t360 * t572 + t366 * t393) * t514 + (t372 * t532 + t373 * t530 + t374 * t528) * t501, m(4) + (t362 * t567 + t368 * t398) * t512 + (t361 * t568 + t367 * t397) * t513 + (t360 * t569 + t366 * t396) * t514 + (t372 * t531 + t373 * t529 + t374 * t527) * t501, t356; t355, t356, 0.3e1 * m(1) + 0.3e1 * m(2) + 0.3e1 * m(3) + m(4);];
MX  = t1;
