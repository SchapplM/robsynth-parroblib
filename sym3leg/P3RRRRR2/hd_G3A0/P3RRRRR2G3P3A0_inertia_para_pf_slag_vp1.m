% Calculate inertia matrix for parallel robot
% P3RRRRR2G3P3A0
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
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:13
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRRRR2G3P3A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR2G3P3A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR2G3P3A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3RRRRR2G3P3A0_inertia_para_pf_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR2G3P3A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRRRR2G3P3A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRRRR2G3P3A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR2G3P3A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR2G3P3A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:11:44
% EndTime: 2020-03-09 21:11:46
% DurationCPUTime: 1.31s
% Computational Cost: add. (2937->207), mult. (6564->416), div. (792->11), fcn. (2934->33), ass. (0->215)
t565 = cos(qJ(3,1));
t682 = t565 ^ 2;
t562 = cos(qJ(3,2));
t681 = t562 ^ 2;
t559 = cos(qJ(3,3));
t680 = t559 ^ 2;
t675 = m(3) / 0.2e1;
t679 = Icges(3,2) / 0.2e1;
t674 = m(3) * rSges(3,2);
t524 = -rSges(3,1) * t674 + Icges(3,4);
t574 = 0.2e1 * qJ(3,3);
t577 = rSges(3,2) ^ 2;
t578 = rSges(3,1) ^ 2;
t672 = (-t577 + t578) * t675 - Icges(3,1) / 0.2e1 + t679;
t678 = cos(t574) * t672 + t524 * sin(t574);
t575 = 0.2e1 * qJ(3,2);
t677 = cos(t575) * t672 + t524 * sin(t575);
t576 = 0.2e1 * qJ(3,1);
t676 = cos(t576) * t672 + t524 * sin(t576);
t673 = m(3) * rSges(3,3);
t561 = cos(qJ(1,3));
t671 = pkin(1) * t561;
t564 = cos(qJ(1,2));
t670 = pkin(1) * t564;
t567 = cos(qJ(1,1));
t669 = pkin(1) * t567;
t547 = legFrame(3,2);
t528 = sin(t547);
t531 = cos(t547);
t550 = sin(qJ(3,3));
t649 = t531 * t550;
t551 = sin(qJ(2,3));
t552 = sin(qJ(1,3));
t560 = cos(qJ(2,3));
t509 = t552 * t551 - t561 * t560;
t662 = t509 * t559;
t488 = t528 * t662 + t649;
t539 = 0.1e1 / t559;
t668 = t488 * t539;
t655 = t528 * t550;
t489 = -t531 * t662 + t655;
t667 = t489 * t539;
t548 = legFrame(2,2);
t529 = sin(t548);
t532 = cos(t548);
t553 = sin(qJ(3,2));
t647 = t532 * t553;
t554 = sin(qJ(2,2));
t555 = sin(qJ(1,2));
t563 = cos(qJ(2,2));
t510 = t555 * t554 - t564 * t563;
t661 = t510 * t562;
t490 = t529 * t661 + t647;
t542 = 0.1e1 / t562;
t666 = t490 * t542;
t653 = t529 * t553;
t491 = -t532 * t661 + t653;
t665 = t491 * t542;
t549 = legFrame(1,2);
t530 = sin(t549);
t533 = cos(t549);
t556 = sin(qJ(3,1));
t645 = t533 * t556;
t557 = sin(qJ(2,1));
t558 = sin(qJ(1,1));
t566 = cos(qJ(2,1));
t511 = t558 * t557 - t567 * t566;
t660 = t511 * t565;
t492 = t530 * t660 + t645;
t545 = 0.1e1 / t565;
t664 = t492 * t545;
t651 = t530 * t556;
t493 = -t533 * t660 + t651;
t663 = t493 * t545;
t525 = sin(qJ(1,3) + qJ(2,3));
t535 = 0.1e1 / t551;
t659 = t525 * t535;
t526 = sin(qJ(1,2) + qJ(2,2));
t536 = 0.1e1 / t554;
t658 = t526 * t536;
t527 = sin(qJ(1,1) + qJ(2,1));
t537 = 0.1e1 / t557;
t657 = t527 * t537;
t656 = t528 * t539;
t654 = t529 * t542;
t652 = t530 * t545;
t650 = t531 * t539;
t648 = t532 * t542;
t646 = t533 * t545;
t644 = t535 * t539;
t540 = 0.1e1 / t680;
t643 = t535 * t540;
t581 = 0.1e1 / pkin(1);
t642 = t535 * t581;
t641 = t536 * t542;
t543 = 0.1e1 / t681;
t640 = t536 * t543;
t639 = t536 * t581;
t638 = t537 * t545;
t546 = 0.1e1 / t682;
t637 = t537 * t546;
t636 = t537 * t581;
t579 = 0.1e1 / pkin(2);
t635 = t539 * t579;
t634 = t540 * t579;
t633 = t542 * t579;
t632 = t543 * t579;
t631 = t545 * t579;
t630 = t546 * t579;
t523 = rSges(3,1) * t673 - Icges(3,5);
t629 = t577 + t578;
t628 = m(3) * rSges(3,1) * pkin(1);
t627 = pkin(1) * t550 * t560;
t626 = pkin(1) * t553 * t563;
t625 = pkin(1) * t556 * t566;
t624 = pkin(2) * t509 * t680;
t623 = pkin(2) * t510 * t681;
t622 = pkin(2) * t511 * t682;
t587 = Icges(2,3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + t679 + Icges(3,1) / 0.2e1;
t588 = 0.2e1 * rSges(3,3) ^ 2 + t629;
t586 = t588 * t675 + t587;
t482 = t586 + t678;
t518 = m(2) * rSges(2,2) - t673;
t571 = m(2) * rSges(2,1);
t585 = (-(-t571 + (-rSges(3,1) * t559 + rSges(3,2) * t550) * m(3)) * t560 - t518 * t551) * pkin(1);
t473 = t585 + t482;
t621 = t473 * t634;
t483 = t586 + t677;
t584 = (-(-t571 + (-rSges(3,1) * t562 + rSges(3,2) * t553) * m(3)) * t563 - t518 * t554) * pkin(1);
t474 = t584 + t483;
t620 = t474 * t632;
t484 = t586 + t676;
t583 = (-(-t571 + (-rSges(3,1) * t565 + rSges(3,2) * t556) * m(3)) * t566 - t518 * t557) * pkin(1);
t475 = t583 + t484;
t619 = t475 * t630;
t476 = -t528 * t624 + (-pkin(2) * t649 + t528 * t671) * t559 - t531 * t627;
t618 = t476 * t643;
t477 = -t529 * t623 + (-pkin(2) * t647 + t529 * t670) * t562 - t532 * t626;
t617 = t477 * t640;
t478 = -t530 * t622 + (-pkin(2) * t645 + t530 * t669) * t565 - t533 * t625;
t616 = t478 * t637;
t479 = t531 * t624 + (-pkin(2) * t655 - t531 * t671) * t559 - t528 * t627;
t615 = t479 * t643;
t480 = t532 * t623 + (-pkin(2) * t653 - t532 * t670) * t562 - t529 * t626;
t614 = t480 * t640;
t481 = t533 * t622 + (-pkin(2) * t651 - t533 * t669) * t565 - t530 * t625;
t613 = t481 * t637;
t612 = t482 * t634;
t611 = t483 * t632;
t610 = t484 * t630;
t609 = t488 * t644;
t608 = t489 * t644;
t607 = t490 * t641;
t606 = t491 * t641;
t605 = t492 * t638;
t604 = t493 * t638;
t494 = pkin(2) * (t561 * t551 + t552 * t560) * t559 + t552 * pkin(1);
t603 = t494 * t644;
t602 = t494 * t635;
t495 = pkin(2) * (t564 * t554 + t555 * t563) * t562 + t555 * pkin(1);
t601 = t495 * t641;
t600 = t495 * t633;
t496 = pkin(2) * (t567 * t557 + t558 * t566) * t565 + t558 * pkin(1);
t599 = t496 * t638;
t598 = t496 * t631;
t522 = -rSges(3,2) * t673 + Icges(3,6);
t497 = t522 * t559 - t523 * t550;
t597 = t497 * t634;
t498 = t522 * t562 - t523 * t553;
t596 = t498 * t632;
t499 = t522 * t565 - t523 * t556;
t595 = t499 * t630;
t594 = t528 * t635;
t593 = t529 * t633;
t592 = t530 * t631;
t591 = t531 * t635;
t590 = t532 * t633;
t589 = t533 * t631;
t580 = pkin(1) ^ 2;
t582 = Icges(1,3) + (0.2e1 * t580 + t588) * t675 + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + m(2) * t580 + t587;
t517 = t629 * m(3) + Icges(3,3);
t487 = (Icges(3,6) + (-pkin(1) * t557 - rSges(3,3)) * t674) * t565 - t556 * (t557 * t628 + t523);
t486 = (Icges(3,6) + (-pkin(1) * t554 - rSges(3,3)) * t674) * t562 - t553 * (t554 * t628 + t523);
t485 = (Icges(3,6) + (-pkin(1) * t551 - rSges(3,3)) * t674) * t559 - t550 * (t551 * t628 + t523);
t472 = (-t487 * t527 + t499 * t598) * t636;
t471 = (-t486 * t526 + t498 * t600) * t639;
t470 = (-t485 * t525 + t497 * t602) * t642;
t469 = t582 + 0.2e1 * t583 + t676;
t468 = t582 + 0.2e1 * t584 + t677;
t467 = t582 + 0.2e1 * t585 + t678;
t466 = t517 * t592 + (t481 * t595 + t487 * t663) * t636;
t465 = t517 * t589 + (t478 * t595 + t487 * t664) * t636;
t464 = t517 * t593 + (t480 * t596 + t486 * t665) * t639;
t463 = t517 * t590 + (t477 * t596 + t486 * t666) * t639;
t462 = t517 * t594 + (t479 * t597 + t485 * t667) * t642;
t461 = t517 * t591 + (t476 * t597 + t485 * t668) * t642;
t460 = (-t475 * t527 + t484 * t598) * t636;
t459 = (-t474 * t526 + t483 * t600) * t639;
t458 = (-t473 * t525 + t482 * t602) * t642;
t457 = (-t469 * t527 + t475 * t598) * t636;
t456 = (-t468 * t526 + t474 * t600) * t639;
t455 = (-t467 * t525 + t473 * t602) * t642;
t454 = t499 * t592 + (t475 * t663 + t481 * t610) * t636;
t453 = t499 * t589 + (t475 * t664 + t478 * t610) * t636;
t452 = t498 * t593 + (t474 * t665 + t480 * t611) * t639;
t451 = t498 * t590 + (t474 * t666 + t477 * t611) * t639;
t450 = t497 * t594 + (t473 * t667 + t479 * t612) * t642;
t449 = t497 * t591 + (t473 * t668 + t476 * t612) * t642;
t448 = t487 * t592 + (t469 * t663 + t481 * t619) * t636;
t447 = t487 * t589 + (t469 * t664 + t478 * t619) * t636;
t446 = t486 * t593 + (t468 * t665 + t480 * t620) * t639;
t445 = t486 * t590 + (t468 * t666 + t477 * t620) * t639;
t444 = t485 * t594 + (t467 * t667 + t479 * t621) * t642;
t443 = t485 * t591 + (t467 * t668 + t476 * t621) * t642;
t1 = [m(4) + (t444 * t608 + t446 * t606 + t448 * t604) * t581 + (t462 * t656 + t464 * t654 + t466 * t652 + (t450 * t615 + t452 * t614 + t454 * t613) * t581) * t579, (t444 * t609 + t446 * t607 + t448 * t605) * t581 + (t462 * t650 + t464 * t648 + t466 * t646 + (t450 * t618 + t452 * t617 + t454 * t616) * t581) * t579, (-t444 * t659 - t446 * t658 - t448 * t657 + (t450 * t603 + t452 * t601 + t454 * t599) * t579) * t581; (t443 * t608 + t445 * t606 + t447 * t604) * t581 + (t461 * t656 + t463 * t654 + t465 * t652 + (t449 * t615 + t451 * t614 + t453 * t613) * t581) * t579, m(4) + (t443 * t609 + t445 * t607 + t447 * t605) * t581 + (t461 * t650 + t463 * t648 + t465 * t646 + (t449 * t618 + t451 * t617 + t453 * t616) * t581) * t579, (-t443 * t659 - t445 * t658 - t447 * t657 + (t449 * t603 + t451 * t601 + t453 * t599) * t579) * t581; (t455 * t608 + t456 * t606 + t457 * t604) * t581 + (t470 * t656 + t471 * t654 + t472 * t652 + (t458 * t615 + t459 * t614 + t460 * t613) * t581) * t579, (t455 * t609 + t456 * t607 + t457 * t605) * t581 + (t470 * t650 + t471 * t648 + t472 * t646 + (t458 * t618 + t459 * t617 + t460 * t616) * t581) * t579, m(4) + (-t455 * t659 - t456 * t658 - t457 * t657 + (t458 * t603 + t459 * t601 + t460 * t599) * t579) * t581;];
MX  = t1;
