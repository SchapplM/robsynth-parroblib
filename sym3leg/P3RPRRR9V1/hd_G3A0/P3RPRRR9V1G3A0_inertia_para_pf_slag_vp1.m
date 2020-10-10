% Calculate inertia matrix for parallel robot
% P3RPRRR9V1G3A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2020-08-06 18:58
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RPRRR9V1G3A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR9V1G3A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR9V1G3A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR9V1G3A0_inertia_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR9V1G3A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRRR9V1G3A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RPRRR9V1G3A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR9V1G3A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR9V1G3A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:56:22
% EndTime: 2020-08-06 18:56:23
% DurationCPUTime: 1.48s
% Computational Cost: add. (4302->275), mult. (6591->480), div. (396->10), fcn. (3540->40), ass. (0->224)
t619 = 2 * pkin(1);
t692 = Icges(3,2) / 0.2e1;
t691 = 2 * pkin(3);
t690 = 4 * rSges(2,3);
t568 = cos(pkin(7));
t689 = 0.2e1 * t568 ^ 2;
t688 = m(2) / 0.2e1;
t687 = rSges(2,2) * m(2);
t594 = rSges(3,2) ^ 2;
t596 = rSges(3,1) ^ 2;
t686 = (-t594 + t596) * m(3) / 0.2e1 - Icges(3,1) / 0.2e1 + t692;
t674 = pkin(5) + qJ(2,3);
t551 = rSges(3,3) + t674;
t685 = m(3) * t551;
t675 = pkin(5) + qJ(2,2);
t552 = rSges(3,3) + t675;
t684 = m(3) * t552;
t676 = pkin(5) + qJ(2,1);
t553 = rSges(3,3) + t676;
t683 = m(3) * t553;
t578 = cos(qJ(3,3));
t564 = t578 ^ 2;
t682 = t564 * pkin(3);
t580 = cos(qJ(3,2));
t565 = t580 ^ 2;
t681 = t565 * pkin(3);
t582 = cos(qJ(3,1));
t566 = t582 ^ 2;
t680 = t566 * pkin(3);
t679 = t578 * pkin(2);
t678 = t580 * pkin(2);
t677 = t582 * pkin(2);
t579 = cos(qJ(1,3));
t554 = pkin(6) + t674;
t573 = sin(qJ(1,3));
t622 = pkin(1) * t579 + t573 * t554;
t567 = sin(pkin(7));
t572 = sin(qJ(3,3));
t628 = t567 * t572;
t478 = t622 * t628 + (t564 - 0.1e1) * t579 * pkin(3);
t487 = pkin(1) * t572 + (-pkin(3) + t679 + 0.2e1 * t682) * t567;
t592 = -pkin(3) / 0.2e1;
t504 = t682 + t679 / 0.2e1 + t592;
t569 = legFrame(3,2);
t542 = sin(t569);
t545 = cos(t569);
t612 = t579 * t628;
t603 = pkin(2) * t612 + (t612 * t691 - t622) * t578;
t538 = pkin(1) * t567;
t625 = t578 * (-t572 * pkin(3) + t538);
t640 = t542 * t579;
t593 = pkin(2) / 0.2e1;
t649 = (t578 * pkin(3) + t593) * t572;
t460 = (-t504 * t640 + t545 * t649) * t689 + (t545 * t487 + t603 * t542) * t568 + t478 * t542 + t545 * t625;
t500 = 0.1e1 / (t568 * t578 - t628);
t673 = t460 * t500;
t634 = t545 * t579;
t461 = (t504 * t634 + t542 * t649) * t689 + (t542 * t487 - t603 * t545) * t568 - t478 * t545 + t542 * t625;
t672 = t461 * t500;
t581 = cos(qJ(1,2));
t555 = pkin(6) + t675;
t575 = sin(qJ(1,2));
t621 = pkin(1) * t581 + t575 * t555;
t574 = sin(qJ(3,2));
t627 = t567 * t574;
t479 = t621 * t627 + (t565 - 0.1e1) * t581 * pkin(3);
t488 = pkin(1) * t574 + (-pkin(3) + t678 + 0.2e1 * t681) * t567;
t505 = t681 + t678 / 0.2e1 + t592;
t570 = legFrame(2,2);
t543 = sin(t570);
t546 = cos(t570);
t611 = t581 * t627;
t602 = pkin(2) * t611 + (t611 * t691 - t621) * t580;
t624 = t580 * (-t574 * pkin(3) + t538);
t638 = t543 * t581;
t648 = (t580 * pkin(3) + t593) * t574;
t462 = (-t505 * t638 + t546 * t648) * t689 + (t546 * t488 + t602 * t543) * t568 + t479 * t543 + t546 * t624;
t501 = 0.1e1 / (t568 * t580 - t627);
t671 = t462 * t501;
t632 = t546 * t581;
t463 = (t505 * t632 + t543 * t648) * t689 + (t543 * t488 - t602 * t546) * t568 - t479 * t546 + t543 * t624;
t670 = t463 * t501;
t583 = cos(qJ(1,1));
t556 = pkin(6) + t676;
t577 = sin(qJ(1,1));
t620 = pkin(1) * t583 + t577 * t556;
t576 = sin(qJ(3,1));
t626 = t567 * t576;
t480 = t620 * t626 + (t566 - 0.1e1) * t583 * pkin(3);
t489 = pkin(1) * t576 + (-pkin(3) + t677 + 0.2e1 * t680) * t567;
t506 = t680 + t677 / 0.2e1 + t592;
t571 = legFrame(1,2);
t544 = sin(t571);
t547 = cos(t571);
t610 = t583 * t626;
t601 = pkin(2) * t610 + (t610 * t691 - t620) * t582;
t623 = t582 * (-t576 * pkin(3) + t538);
t636 = t544 * t583;
t647 = (t582 * pkin(3) + t593) * t576;
t464 = (-t506 * t636 + t547 * t647) * t689 + (t547 * t489 + t601 * t544) * t568 + t480 * t544 + t547 * t623;
t502 = 0.1e1 / (t568 * t582 - t626);
t669 = t464 * t502;
t630 = t547 * t583;
t465 = (t506 * t630 + t544 * t647) * t689 + (t544 * t489 - t601 * t547) * t568 - t480 * t547 + t544 * t623;
t668 = t465 * t502;
t561 = pkin(7) + qJ(3,3);
t532 = sin(t561);
t535 = cos(t561);
t607 = rSges(3,1) * t535 - rSges(3,2) * t532;
t511 = (m(2) * rSges(2,1) + pkin(2) * m(3)) * t568;
t528 = t567 * t687;
t588 = m(2) + m(3);
t608 = -pkin(1) * t588 - t511 + t528;
t475 = -t607 * m(3) + t608;
t522 = t568 * pkin(2) + pkin(1);
t484 = t554 * t579 + (-pkin(3) * t535 - t522) * t573;
t539 = 0.1e1 / t554;
t472 = (-t475 * t573 + t484 * t588) * t539;
t667 = t472 * t500;
t562 = pkin(7) + qJ(3,2);
t533 = sin(t562);
t536 = cos(t562);
t606 = rSges(3,1) * t536 - rSges(3,2) * t533;
t476 = -t606 * m(3) + t608;
t485 = t555 * t581 + (-pkin(3) * t536 - t522) * t575;
t540 = 0.1e1 / t555;
t473 = (-t476 * t575 + t485 * t588) * t540;
t666 = t473 * t501;
t563 = pkin(7) + qJ(3,1);
t534 = sin(t563);
t537 = cos(t563);
t605 = rSges(3,1) * t537 - rSges(3,2) * t534;
t477 = -t605 * m(3) + t608;
t486 = t556 * t583 + (-pkin(3) * t537 - t522) * t577;
t541 = 0.1e1 / t556;
t474 = (-t477 * t577 + t486 * t588) * t541;
t665 = t474 * t502;
t490 = t545 * t532 - t535 * t640;
t529 = 0.1e1 / t535;
t664 = t490 * t529;
t663 = t490 * t539;
t491 = t542 * t532 + t535 * t634;
t662 = t491 * t529;
t661 = t491 * t539;
t492 = t546 * t533 - t536 * t638;
t530 = 0.1e1 / t536;
t660 = t492 * t530;
t659 = t492 * t540;
t493 = t543 * t533 + t536 * t632;
t658 = t493 * t530;
t657 = t493 * t540;
t494 = t547 * t534 - t537 * t636;
t531 = 0.1e1 / t537;
t656 = t494 * t531;
t655 = t494 * t541;
t495 = t544 * t534 + t537 * t630;
t654 = t495 * t531;
t653 = t495 * t541;
t652 = t500 * t588;
t651 = t501 * t588;
t650 = t502 * t588;
t646 = t529 * t542;
t645 = t529 * t545;
t644 = t530 * t543;
t643 = t530 * t546;
t642 = t531 * t544;
t641 = t531 * t547;
t598 = 1 / pkin(3);
t639 = t542 * t598;
t637 = t543 * t598;
t635 = t544 * t598;
t633 = t545 * t598;
t631 = t546 * t598;
t629 = t547 * t598;
t618 = t475 * t500 * t539;
t617 = t476 * t501 * t540;
t616 = t477 * t502 * t541;
t481 = (-rSges(3,2) * t685 + Icges(3,6)) * t535 - (rSges(3,1) * t685 - Icges(3,5)) * t532;
t615 = t481 * t573 * t598;
t482 = (-rSges(3,2) * t684 + Icges(3,6)) * t536 - (rSges(3,1) * t684 - Icges(3,5)) * t533;
t614 = t482 * t575 * t598;
t483 = (-rSges(3,2) * t683 + Icges(3,6)) * t537 - (rSges(3,1) * t683 - Icges(3,5)) * t534;
t613 = t483 * t577 * t598;
t590 = 2 * pkin(1) ^ 2;
t595 = rSges(2,2) ^ 2;
t597 = rSges(2,1) ^ 2;
t609 = (2 * rSges(2,3) ^ 2) + t590 + t595 + t597;
t599 = pkin(2) ^ 2;
t604 = t590 / 0.2e1 + t594 / 0.2e1 + t596 / 0.2e1 + t599 / 0.2e1;
t591 = 0.2e1 * pkin(7);
t600 = Icges(1,3) + (m(3) * t599 + (-t595 + t597) * m(2) - Icges(2,1) + Icges(2,2)) * cos(t591) / 0.2e1 + (-rSges(2,1) * t687 + Icges(2,4)) * sin(t591) + t511 * t619 + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) - 0.2e1 * pkin(1) * t528 + t692 + Icges(2,2) / 0.2e1 + Icges(3,1) / 0.2e1 + Icges(2,1) / 0.2e1;
t560 = t591 + qJ(3,1);
t559 = t591 + qJ(3,2);
t558 = t591 + qJ(3,3);
t527 = -m(3) * rSges(3,1) * rSges(3,2) + Icges(3,4);
t525 = 0.2e1 * t563;
t524 = 0.2e1 * t562;
t523 = 0.2e1 * t561;
t514 = (t594 + t596) * m(3) + Icges(3,3);
t471 = (t483 * t653 + t514 * t635) * t531;
t470 = (t483 * t655 + t514 * t629) * t531;
t469 = (t482 * t657 + t514 * t637) * t530;
t468 = (t482 * t659 + t514 * t631) * t530;
t467 = (t481 * t661 + t514 * t639) * t529;
t466 = (t481 * t663 + t514 * t633) * t529;
t459 = t527 * sin(t525) + (t553 ^ 2 + t605 * t619 + ((-sin(t560) - t576) * rSges(3,2) + (cos(t560) + t582) * rSges(3,1)) * pkin(2) + t604) * m(3) + (((t690 + 2 * qJ(2,1)) * qJ(2,1)) + t609) * t688 + t600 + cos(t525) * t686;
t458 = t527 * sin(t524) + (t552 ^ 2 + t606 * t619 + ((-sin(t559) - t574) * rSges(3,2) + (cos(t559) + t580) * rSges(3,1)) * pkin(2) + t604) * m(3) + (((t690 + 2 * qJ(2,2)) * qJ(2,2)) + t609) * t688 + t600 + cos(t524) * t686;
t457 = t527 * sin(t523) + (t551 ^ 2 + t607 * t619 + ((-sin(t558) - t572) * rSges(3,2) + (cos(t558) + t578) * rSges(3,1)) * pkin(2) + t604) * m(3) + (((t690 + 2 * qJ(2,3)) * qJ(2,3)) + t609) * t688 + t600 + cos(t523) * t686;
t456 = (-t459 * t577 + t477 * t486) * t541;
t455 = (-t458 * t575 + t476 * t485) * t540;
t454 = (-t457 * t573 + t475 * t484) * t539;
t453 = (t465 * t650 + t477 * t654) * t541;
t452 = (t464 * t650 + t477 * t656) * t541;
t451 = (t463 * t651 + t476 * t658) * t540;
t450 = (t462 * t651 + t476 * t660) * t540;
t449 = (t461 * t652 + t475 * t662) * t539;
t448 = (t460 * t652 + t475 * t664) * t539;
t447 = t465 * t616 + (t459 * t653 + t483 * t635) * t531;
t446 = t464 * t616 + (t459 * t655 + t483 * t629) * t531;
t445 = t463 * t617 + (t458 * t657 + t482 * t637) * t530;
t444 = t462 * t617 + (t458 * t659 + t482 * t631) * t530;
t443 = t461 * t618 + (t457 * t661 + t481 * t639) * t529;
t442 = t460 * t618 + (t457 * t663 + t481 * t633) * t529;
t1 = [m(4) + (t447 * t654 + t453 * t668) * t541 + (t445 * t658 + t451 * t670) * t540 + (t443 * t662 + t449 * t672) * t539 + (t467 * t646 + t469 * t644 + t471 * t642) * t598, (t447 * t656 + t453 * t669) * t541 + (t445 * t660 + t451 * t671) * t540 + (t443 * t664 + t449 * t673) * t539 + (t467 * t645 + t469 * t643 + t471 * t641) * t598, (-t447 * t577 + t453 * t486) * t541 + (-t445 * t575 + t451 * t485) * t540 + (-t443 * t573 + t449 * t484) * t539; (t446 * t654 + t452 * t668) * t541 + (t444 * t658 + t450 * t670) * t540 + (t442 * t662 + t448 * t672) * t539 + (t466 * t646 + t468 * t644 + t470 * t642) * t598, m(4) + (t446 * t656 + t452 * t669) * t541 + (t444 * t660 + t450 * t671) * t540 + (t442 * t664 + t448 * t673) * t539 + (t466 * t645 + t468 * t643 + t470 * t641) * t598, (-t446 * t577 + t452 * t486) * t541 + (-t444 * t575 + t450 * t485) * t540 + (-t442 * t573 + t448 * t484) * t539; (t465 * t665 + (t456 * t495 - t544 * t613) * t531) * t541 + (t463 * t666 + (t455 * t493 - t543 * t614) * t530) * t540 + (t461 * t667 + (t454 * t491 - t542 * t615) * t529) * t539, (t464 * t665 + (t456 * t494 - t547 * t613) * t531) * t541 + (t462 * t666 + (t455 * t492 - t546 * t614) * t530) * t540 + (t460 * t667 + (t454 * t490 - t545 * t615) * t529) * t539, m(4) + (-t456 * t577 + t474 * t486) * t541 + (-t455 * t575 + t473 * t485) * t540 + (-t454 * t573 + t472 * t484) * t539;];
MX  = t1;
