% Calculate inertia matrix for parallel robot
% P3PRRRR8V1G2A0
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
% Datum: 2020-08-06 17:03
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3PRRRR8V1G2A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G2A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G2A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G2A0_inertia_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V1G2A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRRR8V1G2A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRRR8V1G2A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G2A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G2A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:02:35
% EndTime: 2020-08-06 17:02:37
% DurationCPUTime: 1.42s
% Computational Cost: add. (3072->211), mult. (8334->417), div. (648->7), fcn. (7470->28), ass. (0->206)
t610 = m(3) / 0.2e1;
t609 = Icges(3,2) / 0.2e1;
t492 = sin(qJ(2,3));
t498 = cos(qJ(2,3));
t497 = cos(qJ(3,3));
t601 = pkin(2) * t497;
t457 = -t498 * pkin(5) + t492 * t601;
t485 = sin(pkin(3));
t487 = cos(pkin(3));
t491 = sin(qJ(3,3));
t570 = t491 * t487;
t451 = pkin(2) * t570 + t457 * t485;
t608 = 0.1e1 / t451;
t494 = sin(qJ(2,2));
t500 = cos(qJ(2,2));
t499 = cos(qJ(3,2));
t600 = pkin(2) * t499;
t458 = -t500 * pkin(5) + t494 * t600;
t493 = sin(qJ(3,2));
t568 = t493 * t487;
t452 = pkin(2) * t568 + t458 * t485;
t607 = 0.1e1 / t452;
t496 = sin(qJ(2,1));
t502 = cos(qJ(2,1));
t501 = cos(qJ(3,1));
t599 = pkin(2) * t501;
t459 = -t502 * pkin(5) + t496 * t599;
t495 = sin(qJ(3,1));
t566 = t495 * t487;
t453 = pkin(2) * t566 + t459 * t485;
t606 = 0.1e1 / t453;
t605 = m(3) * rSges(3,3);
t509 = rSges(3,2) ^ 2;
t510 = rSges(3,1) ^ 2;
t604 = (-t509 + t510) * t610 - Icges(3,1) / 0.2e1 + t609;
t511 = 0.1e1 / pkin(2);
t603 = m(3) * t511;
t602 = pkin(2) * t485;
t460 = pkin(5) * t492 + t498 * t601;
t484 = sin(pkin(6));
t486 = cos(pkin(6));
t526 = -t457 * t487 + t491 * t602;
t421 = -t484 * t460 + t526 * t486;
t488 = legFrame(3,2);
t474 = sin(t488);
t477 = cos(t488);
t412 = t421 * t474 + t451 * t477;
t598 = t412 * t608;
t413 = -t421 * t477 + t451 * t474;
t597 = t413 * t608;
t461 = pkin(5) * t494 + t500 * t600;
t525 = -t458 * t487 + t493 * t602;
t422 = -t484 * t461 + t525 * t486;
t489 = legFrame(2,2);
t475 = sin(t489);
t478 = cos(t489);
t414 = t422 * t475 + t452 * t478;
t596 = t414 * t607;
t415 = -t422 * t478 + t452 * t475;
t595 = t415 * t607;
t462 = pkin(5) * t496 + t502 * t599;
t524 = -t459 * t487 + t495 * t602;
t423 = -t484 * t462 + t524 * t486;
t490 = legFrame(1,2);
t476 = sin(t490);
t479 = cos(t490);
t416 = t423 * t476 + t453 * t479;
t594 = t416 * t606;
t417 = -t423 * t479 + t453 * t476;
t593 = t417 * t606;
t424 = t486 * t460 + t526 * t484;
t592 = t424 * t608;
t425 = t486 * t461 + t525 * t484;
t591 = t425 * t607;
t426 = t486 * t462 + t524 * t484;
t590 = t426 * t606;
t470 = m(2) * rSges(2,2) - t605;
t505 = m(2) * rSges(2,1);
t530 = rSges(3,1) * t497 - rSges(3,2) * t491;
t589 = ((t530 * m(3) + t505) * t498 - t470 * t492) * t485;
t529 = rSges(3,1) * t499 - rSges(3,2) * t493;
t588 = ((t529 * m(3) + t505) * t500 - t470 * t494) * t485;
t528 = rSges(3,1) * t501 - rSges(3,2) * t495;
t587 = ((t528 * m(3) + t505) * t502 - t470 * t496) * t485;
t586 = t608 / t497;
t585 = t607 / t499;
t584 = t606 / t501;
t480 = m(1) + m(2) + m(3);
t583 = t608 * t480;
t582 = t607 * t480;
t581 = t606 * t480;
t471 = -rSges(3,2) * t605 + Icges(3,6);
t472 = rSges(3,1) * t605 - Icges(3,5);
t454 = t471 * t497 - t472 * t491;
t580 = t454 * t511;
t455 = t471 * t499 - t472 * t493;
t579 = t455 * t511;
t456 = t471 * t501 - t472 * t495;
t578 = t456 * t511;
t564 = t509 + t510;
t577 = (t564 * m(3) + Icges(3,3)) * t511;
t576 = t487 * t492;
t575 = t487 * t494;
t574 = t487 * t496;
t573 = t487 * t498;
t572 = t487 * t500;
t571 = t487 * t502;
t569 = t491 * t498;
t567 = t493 * t500;
t565 = t495 * t502;
t439 = t530 * t487 - t485 * t492 * (t491 * rSges(3,1) + t497 * rSges(3,2));
t563 = m(3) * t439 * t608;
t562 = t439 * t603;
t440 = t529 * t487 - t485 * t494 * (t493 * rSges(3,1) + t499 * rSges(3,2));
t561 = m(3) * t440 * t607;
t560 = t440 * t603;
t441 = t528 * t487 - t485 * t496 * (t495 * rSges(3,1) + t501 * rSges(3,2));
t559 = m(3) * t441 * t606;
t558 = t441 * t603;
t427 = -(-t484 * t492 + t486 * t573) * t601 - pkin(5) * (t484 * t498 + t486 * t576);
t557 = t427 * t586;
t428 = -(-t484 * t494 + t486 * t572) * t600 - pkin(5) * (t484 * t500 + t486 * t575);
t556 = t428 * t585;
t429 = -(-t484 * t496 + t486 * t571) * t599 - pkin(5) * (t484 * t502 + t486 * t574);
t555 = t429 * t584;
t517 = t485 * t497 + t492 * t570;
t436 = -t484 * t569 - t517 * t486;
t554 = t436 * t586;
t516 = t485 * t499 + t494 * t568;
t437 = -t484 * t567 - t516 * t486;
t553 = t437 * t585;
t515 = t485 * t501 + t496 * t566;
t438 = -t484 * t565 - t515 * t486;
t552 = t438 * t584;
t551 = t608 * t589;
t550 = t607 * t588;
t549 = t606 * t587;
t548 = t474 * t586;
t547 = t477 * t586;
t546 = t475 * t585;
t545 = t478 * t585;
t544 = t476 * t584;
t543 = t479 * t584;
t430 = (t484 * t573 + t486 * t492) * t601 + (t484 * t576 - t486 * t498) * pkin(5);
t542 = t430 * t548;
t541 = t430 * t547;
t431 = (t484 * t572 + t486 * t494) * t600 + (t484 * t575 - t486 * t500) * pkin(5);
t540 = t431 * t546;
t539 = t431 * t545;
t432 = (t484 * t571 + t486 * t496) * t599 + (t484 * t574 - t486 * t502) * pkin(5);
t538 = t432 * t544;
t537 = t432 * t543;
t433 = t517 * t484 - t486 * t569;
t536 = t433 * t548;
t535 = t433 * t547;
t434 = t516 * t484 - t486 * t567;
t534 = t434 * t546;
t533 = t434 * t545;
t435 = t515 * t484 - t486 * t565;
t532 = t435 * t544;
t531 = t435 * t543;
t527 = Icges(2,3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + (0.2e1 * rSges(3,3) ^ 2 + t564) * t610 + t609 + Icges(3,1) / 0.2e1;
t473 = -m(3) * rSges(3,1) * rSges(3,2) + Icges(3,4);
t506 = 0.2e1 * qJ(3,3);
t418 = cos(t506) * t604 + t473 * sin(t506) + t527;
t523 = t418 * t433 + t430 * t580;
t507 = 0.2e1 * qJ(3,2);
t419 = cos(t507) * t604 + t473 * sin(t507) + t527;
t522 = t419 * t434 + t431 * t579;
t508 = 0.2e1 * qJ(3,1);
t420 = cos(t508) * t604 + t473 * sin(t508) + t527;
t521 = t420 * t435 + t432 * t578;
t520 = t430 * t577 + t433 * t454;
t519 = t431 * t577 + t434 * t455;
t518 = t432 * t577 + t435 * t456;
t514 = t430 * t562 + t433 * t589;
t513 = t431 * t560 + t434 * t588;
t512 = t432 * t558 + t435 * t587;
t411 = t426 * t559 + (t429 * t577 + t438 * t456) * t584;
t410 = t425 * t561 + (t428 * t577 + t437 * t455) * t585;
t409 = t424 * t563 + (t427 * t577 + t436 * t454) * t586;
t408 = t426 * t581 + (t429 * t558 + t438 * t587) * t584;
t407 = t425 * t582 + (t428 * t560 + t437 * t588) * t585;
t406 = t424 * t583 + (t427 * t562 + t436 * t589) * t586;
t405 = t426 * t549 + (t420 * t438 + t429 * t578) * t584;
t404 = t425 * t550 + (t419 * t437 + t428 * t579) * t585;
t403 = t424 * t551 + (t418 * t436 + t427 * t580) * t586;
t402 = t417 * t559 - t518 * t543;
t401 = t416 * t559 + t518 * t544;
t400 = t415 * t561 - t519 * t545;
t399 = t414 * t561 + t519 * t546;
t398 = t413 * t563 - t520 * t547;
t397 = t412 * t563 + t520 * t548;
t396 = t417 * t581 - t512 * t543;
t395 = t416 * t581 + t512 * t544;
t394 = t415 * t582 - t513 * t545;
t393 = t414 * t582 + t513 * t546;
t392 = t413 * t583 - t514 * t547;
t391 = t412 * t583 + t514 * t548;
t390 = t417 * t549 - t521 * t543;
t389 = t416 * t549 + t521 * t544;
t388 = t415 * t550 - t522 * t545;
t387 = t414 * t550 + t522 * t546;
t386 = t413 * t551 - t523 * t547;
t385 = t412 * t551 + t523 * t548;
t1 = [-t386 * t535 - t388 * t533 - t390 * t531 + t392 * t597 + t394 * t595 + t396 * t593 + m(4) + (-t398 * t541 - t400 * t539 - t402 * t537) * t511, t386 * t536 + t388 * t534 + t390 * t532 + t392 * t598 + t394 * t596 + t396 * t594 + (t398 * t542 + t400 * t540 + t402 * t538) * t511, t386 * t554 + t388 * t553 + t390 * t552 + t392 * t592 + t394 * t591 + t396 * t590 + (t398 * t557 + t400 * t556 + t402 * t555) * t511; -t385 * t535 - t387 * t533 - t389 * t531 + t391 * t597 + t393 * t595 + t395 * t593 + (-t397 * t541 - t399 * t539 - t401 * t537) * t511, t385 * t536 + t387 * t534 + t389 * t532 + t391 * t598 + t393 * t596 + t395 * t594 + m(4) + (t397 * t542 + t399 * t540 + t401 * t538) * t511, t385 * t554 + t387 * t553 + t389 * t552 + t391 * t592 + t393 * t591 + t395 * t590 + (t397 * t557 + t399 * t556 + t401 * t555) * t511; -t403 * t535 - t404 * t533 - t405 * t531 + t406 * t597 + t407 * t595 + t408 * t593 + (-t409 * t541 - t410 * t539 - t411 * t537) * t511, t403 * t536 + t404 * t534 + t405 * t532 + t406 * t598 + t407 * t596 + t408 * t594 + (t409 * t542 + t410 * t540 + t411 * t538) * t511, t403 * t554 + t404 * t553 + t405 * t552 + t406 * t592 + t407 * t591 + t408 * t590 + m(4) + (t409 * t557 + t410 * t556 + t411 * t555) * t511;];
MX  = t1;
