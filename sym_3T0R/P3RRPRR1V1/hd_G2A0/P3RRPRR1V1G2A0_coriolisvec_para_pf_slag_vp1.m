% Calculate vector of centrifugal and Coriolis load for parallel robot
% P3RRPRR1V1G2A0
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
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4]';
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
% Datum: 2022-11-04 17:08
% Revision: e482436b586c4f286726c907c195760c5ac72455 (2022-11-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RRPRR1V1G2A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR1V1G2A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR1V1G2A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR1V1G2A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3RRPRR1V1G2A0_coriolisvec_para_pf_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR1V1G2A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR1V1G2A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRPRR1V1G2A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR1V1G2A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR1V1G2A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-11-04 17:07:59
% EndTime: 2022-11-04 17:08:00
% DurationCPUTime: 1.00s
% Computational Cost: add. (5685->185), mult. (6306->365), div. (2331->7), fcn. (5724->24), ass. (0->152)
t525 = legFrame(3,2);
t498 = sin(t525);
t501 = cos(t525);
t534 = cos(qJ(2,3));
t512 = 0.1e1 / t534;
t548 = pkin(2) + pkin(1);
t518 = 0.1e1 / t548;
t546 = xDP(2);
t547 = xDP(1);
t467 = (t498 * t547 + t501 * t546) * t518 * t512;
t464 = t467 ^ 2;
t526 = legFrame(2,2);
t499 = sin(t526);
t502 = cos(t526);
t536 = cos(qJ(2,2));
t514 = 0.1e1 / t536;
t468 = (t499 * t547 + t502 * t546) * t518 * t514;
t465 = t468 ^ 2;
t527 = legFrame(1,2);
t500 = sin(t527);
t503 = cos(t527);
t538 = cos(qJ(2,1));
t516 = 0.1e1 / t538;
t469 = (t500 * t547 + t503 * t546) * t518 * t516;
t466 = t469 ^ 2;
t517 = t548 ^ 2;
t528 = sin(qJ(2,3));
t529 = sin(qJ(1,3));
t592 = t529 * t534;
t476 = -t498 * t592 + t528 * t501;
t479 = t498 * t528 + t501 * t592;
t522 = pkin(3) + qJ(3,3);
t507 = 0.1e1 / t522;
t535 = cos(qJ(1,3));
t545 = xDP(3);
t446 = (t535 * t545 + (t476 * t546 + t479 * t547) * t512) * t507;
t612 = 0.2e1 * t446;
t530 = sin(qJ(2,2));
t531 = sin(qJ(1,2));
t590 = t531 * t536;
t477 = -t499 * t590 + t530 * t502;
t480 = t499 * t530 + t502 * t590;
t523 = pkin(3) + qJ(3,2);
t508 = 0.1e1 / t523;
t537 = cos(qJ(1,2));
t447 = (t537 * t545 + (t477 * t546 + t480 * t547) * t514) * t508;
t611 = 0.2e1 * t447;
t532 = sin(qJ(2,1));
t533 = sin(qJ(1,1));
t588 = t533 * t538;
t478 = -t500 * t588 + t532 * t503;
t481 = t500 * t532 + t503 * t588;
t524 = pkin(3) + qJ(3,1);
t509 = 0.1e1 / t524;
t539 = cos(qJ(1,1));
t448 = (t539 * t545 + (t478 * t546 + t481 * t547) * t516) * t509;
t610 = 0.2e1 * t448;
t609 = -m(3) / 0.2e1;
t608 = rSges(2,3) * m(2);
t540 = pkin(1) + rSges(3,1);
t553 = rSges(2,2) ^ 2;
t555 = rSges(2,1) ^ 2;
t557 = (rSges(3,2) + t540) * (-rSges(3,2) + t540) * m(3) - Icges(2,1) - Icges(3,1) + Icges(2,2) + Icges(3,2) + (-t553 + t555) * m(2);
t607 = -t557 / 0.2e1;
t519 = rSges(3,3) + qJ(3,3);
t606 = m(3) * t519;
t520 = rSges(3,3) + qJ(3,2);
t605 = m(3) * t520;
t521 = rSges(3,3) + qJ(3,1);
t604 = m(3) * t521;
t603 = m(3) * t540;
t602 = t464 * t528;
t601 = t465 * t530;
t600 = t466 * t532;
t599 = t557 * t528;
t598 = t557 * t530;
t597 = t557 * t532;
t491 = m(2) * rSges(2,1) * rSges(2,2) + rSges(3,2) * t603 - Icges(2,4) - Icges(3,4);
t596 = t491 * t467;
t595 = t491 * t468;
t594 = t491 * t469;
t593 = t528 * t548;
t591 = t530 * t548;
t589 = t532 * t548;
t587 = t548 * t534;
t586 = t548 * t536;
t585 = t548 * t538;
t561 = -t522 * t535 + t529 * t587;
t458 = -t561 * t498 + t501 * t593;
t461 = t498 * t593 + t561 * t501;
t482 = t529 * t522 + t535 * t587;
t443 = (t458 * t546 + t461 * t547 + t482 * t545) * t507;
t511 = t534 ^ 2;
t580 = 0.2e1 * t548;
t425 = (-t517 * t464 + ((-t517 * t511 - t522 ^ 2) * t446 + (t467 * t528 * t522 + t443 * t534) * t580) * t446) * t507;
t437 = (-t548 * t464 * t512 + (-t446 * t587 + 0.2e1 * t443) * t446) * t507;
t565 = rSges(2,1) * t608 - Icges(2,5) - Icges(3,5);
t485 = t519 * t603 + t565;
t576 = rSges(2,2) * t608 - Icges(2,6) - Icges(3,6);
t488 = rSges(3,2) * t606 + t576;
t455 = t485 * t528 + t534 * t488;
t492 = -t528 * rSges(3,2) + t540 * t534;
t549 = 0.2e1 * qJ(2,3);
t581 = t553 + t555;
t558 = -(rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) - ((2 * rSges(2,3) ^ 2) + t581) * m(2) / 0.2e1 - Icges(3,2) / 0.2e1 - Icges(2,2) / 0.2e1 - Icges(3,1) / 0.2e1 - Icges(2,1) / 0.2e1 - Icges(1,3);
t572 = rSges(3,2) ^ 2 + (pkin(1) ^ 2) + ((2 * pkin(1) + rSges(3,1)) * rSges(3,1));
t579 = t512 * t602;
t584 = (cos(t549) * t607 + t491 * sin(t549) + (0.2e1 * t519 ^ 2 + t572) * t609 + t558) * t437 - t455 * t579 + t492 * m(3) * t425 - 0.4e1 * t446 * t511 * t596 - 0.2e1 * t467 * (t446 * t599 + t485 * t467 / 0.2e1) * t534 + t488 * t602 + (t443 * t606 + t596) * t612;
t560 = -t523 * t537 + t531 * t586;
t459 = -t560 * t499 + t502 * t591;
t462 = t499 * t591 + t560 * t502;
t483 = t531 * t523 + t537 * t586;
t444 = (t459 * t546 + t462 * t547 + t483 * t545) * t508;
t513 = t536 ^ 2;
t426 = (-t517 * t465 + ((-t517 * t513 - t523 ^ 2) * t447 + (t468 * t530 * t523 + t444 * t536) * t580) * t447) * t508;
t438 = (-t548 * t465 * t514 + (-t447 * t586 + 0.2e1 * t444) * t447) * t508;
t486 = t520 * t603 + t565;
t489 = rSges(3,2) * t605 + t576;
t456 = t486 * t530 + t536 * t489;
t493 = -t530 * rSges(3,2) + t540 * t536;
t550 = 0.2e1 * qJ(2,2);
t578 = t514 * t601;
t583 = (cos(t550) * t607 + t491 * sin(t550) + (0.2e1 * t520 ^ 2 + t572) * t609 + t558) * t438 - t456 * t578 + t493 * m(3) * t426 - 0.4e1 * t447 * t513 * t595 - 0.2e1 * t468 * (t447 * t598 + t486 * t468 / 0.2e1) * t536 + t489 * t601 + (t444 * t605 + t595) * t611;
t559 = -t524 * t539 + t533 * t585;
t460 = -t559 * t500 + t503 * t589;
t463 = t500 * t589 + t559 * t503;
t484 = t533 * t524 + t539 * t585;
t445 = (t460 * t546 + t463 * t547 + t484 * t545) * t509;
t515 = t538 ^ 2;
t427 = (-t517 * t466 + ((-t517 * t515 - t524 ^ 2) * t448 + (t469 * t532 * t524 + t445 * t538) * t580) * t448) * t509;
t439 = (-t548 * t466 * t516 + (-t448 * t585 + 0.2e1 * t445) * t448) * t509;
t487 = t521 * t603 + t565;
t490 = rSges(3,2) * t604 + t576;
t457 = t487 * t532 + t538 * t490;
t494 = -t532 * rSges(3,2) + t540 * t538;
t551 = 0.2e1 * qJ(2,1);
t577 = t516 * t600;
t582 = (cos(t551) * t607 + t491 * sin(t551) + (0.2e1 * t521 ^ 2 + t572) * t609 + t558) * t439 - t457 * t577 + t494 * m(3) * t427 - 0.4e1 * t448 * t515 * t594 - 0.2e1 * t469 * (t448 * t597 + t487 * t469 / 0.2e1) * t538 + t490 * t600 + (t445 * t604 + t594) * t610;
t575 = t584 * t512;
t574 = t583 * t514;
t573 = t582 * t516;
t568 = rSges(3,2) * t534 + t540 * t528;
t571 = ((-t446 * t519 / 0.2e1 + t568 * t467) * t612 + t437 * t492 - t425) * m(3);
t567 = rSges(3,2) * t536 + t540 * t530;
t570 = ((-t447 * t520 / 0.2e1 + t567 * t468) * t611 + t438 * t493 - t426) * m(3);
t566 = rSges(3,2) * t538 + t540 * t532;
t569 = ((-t448 * t521 / 0.2e1 + t566 * t469) * t610 + t439 * t494 - t427) * m(3);
t475 = -t581 * m(2) - t572 * m(3) - Icges(2,3) - Icges(3,3);
t564 = ((-t568 * t443 * m(3) + (t534 * t599 / 0.2e1 + (t511 - 0.1e1 / 0.2e1) * t491) * t446) * t612 + t455 * t437 - t475 * t579) * t512;
t563 = ((-t567 * t444 * m(3) + (t536 * t598 / 0.2e1 + (t513 - 0.1e1 / 0.2e1) * t491) * t447) * t611 + t456 * t438 - t475 * t578) * t514;
t562 = ((-t566 * t445 * m(3) + (t538 * t597 / 0.2e1 + (t515 - 0.1e1 / 0.2e1) * t491) * t448) * t610 + t457 * t439 - t475 * t577) * t516;
t1 = [(t569 * t463 + t481 * t573) * t509 + (t570 * t462 + t480 * t574) * t508 + (t571 * t461 + t479 * t575) * t507 + (t498 * t564 + t499 * t563 + t500 * t562) * t518; (t569 * t460 + t478 * t573) * t509 + (t570 * t459 + t477 * t574) * t508 + (t571 * t458 + t476 * t575) * t507 + (t501 * t564 + t502 * t563 + t503 * t562) * t518; (t569 * t484 + t582 * t539) * t509 + (t570 * t483 + t583 * t537) * t508 + (t571 * t482 + t584 * t535) * t507;];
taucX  = t1;
