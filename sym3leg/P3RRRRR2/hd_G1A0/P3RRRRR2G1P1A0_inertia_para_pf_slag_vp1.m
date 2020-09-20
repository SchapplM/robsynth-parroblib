% Calculate inertia matrix for parallel robot
% P3RRRRR2G1P1A0
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
% Datum: 2020-03-09 21:05
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRRRR2G1P1A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR2G1P1A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR2G1P1A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3RRRRR2G1P1A0_inertia_para_pf_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR2G1P1A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRRRR2G1P1A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRRRR2G1P1A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR2G1P1A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR2G1P1A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:04:07
% EndTime: 2020-03-09 21:04:08
% DurationCPUTime: 0.96s
% Computational Cost: add. (3105->189), mult. (5004->311), div. (570->14), fcn. (1338->48), ass. (0->164)
t643 = m(3) / 0.2e1;
t651 = Icges(3,2) / 0.2e1;
t642 = m(3) * rSges(3,2);
t537 = -rSges(3,1) * t642 + Icges(3,4);
t572 = 2 * qJ(3,3);
t575 = rSges(3,2) ^ 2;
t576 = rSges(3,1) ^ 2;
t640 = (-t575 + t576) * t643 - Icges(3,1) / 0.2e1 + t651;
t650 = cos(t572) * t640 + t537 * sin(t572);
t573 = 2 * qJ(3,2);
t649 = cos(t573) * t640 + t537 * sin(t573);
t574 = 2 * qJ(3,1);
t648 = cos(t574) * t640 + t537 * sin(t574);
t647 = -2 * pkin(1);
t560 = cos(qJ(3,3));
t545 = 0.1e1 / t560;
t646 = t545 ^ 2;
t562 = cos(qJ(3,2));
t547 = 0.1e1 / t562;
t645 = t547 ^ 2;
t564 = cos(qJ(3,1));
t549 = 0.1e1 / t564;
t644 = t549 ^ 2;
t641 = m(3) * rSges(3,3);
t536 = rSges(3,1) * t641 - Icges(3,5);
t554 = sin(qJ(3,3));
t555 = sin(qJ(2,3));
t611 = m(3) * rSges(3,1) * pkin(1);
t486 = (Icges(3,6) + (-pkin(1) * t555 - rSges(3,3)) * t642) * t560 - t554 * (t555 * t611 + t536);
t639 = t486 * t545;
t556 = sin(qJ(3,2));
t557 = sin(qJ(2,2));
t487 = (Icges(3,6) + (-pkin(1) * t557 - rSges(3,3)) * t642) * t562 - t556 * (t557 * t611 + t536);
t638 = t487 * t547;
t558 = sin(qJ(3,1));
t559 = sin(qJ(2,1));
t488 = (Icges(3,6) + (-pkin(1) * t559 - rSges(3,3)) * t642) * t564 - t558 * (t559 * t611 + t536);
t637 = t488 * t549;
t538 = qJ(1,3) + legFrame(3,3);
t531 = qJ(2,3) + t538;
t525 = qJ(3,3) + t531;
t526 = -qJ(3,3) + t531;
t489 = sin(t538) * t647 + (-sin(t526) - sin(t525)) * pkin(2);
t508 = 0.1e1 / (sin(qJ(2,3) + qJ(3,3)) + sin(qJ(2,3) - qJ(3,3)));
t636 = t489 * t508;
t539 = qJ(1,2) + legFrame(2,3);
t532 = qJ(2,2) + t539;
t527 = qJ(3,2) + t532;
t528 = -qJ(3,2) + t532;
t490 = sin(t539) * t647 + (-sin(t528) - sin(t527)) * pkin(2);
t509 = 0.1e1 / (sin(qJ(2,2) + qJ(3,2)) + sin(qJ(2,2) - qJ(3,2)));
t635 = t490 * t509;
t540 = qJ(1,1) + legFrame(1,3);
t533 = qJ(2,1) + t540;
t529 = qJ(3,1) + t533;
t530 = -qJ(3,1) + t533;
t491 = sin(t540) * t647 + (-sin(t530) - sin(t529)) * pkin(2);
t510 = 0.1e1 / (sin(qJ(2,1) + qJ(3,1)) + sin(qJ(2,1) - qJ(3,1)));
t634 = t491 * t510;
t492 = cos(t538) * t647 + (-cos(t526) - cos(t525)) * pkin(2);
t633 = t492 * t508;
t493 = cos(t539) * t647 + (-cos(t528) - cos(t527)) * pkin(2);
t632 = t493 * t509;
t494 = cos(t540) * t647 + (-cos(t530) - cos(t529)) * pkin(2);
t631 = t494 * t510;
t577 = 0.1e1 / pkin(2);
t630 = t508 * t577;
t629 = t509 * t577;
t628 = t510 * t577;
t561 = cos(qJ(2,3));
t627 = (pkin(1) * t561 + pkin(2) * t560) / t560 ^ 2;
t563 = cos(qJ(2,2));
t626 = (pkin(1) * t563 + pkin(2) * t562) / t562 ^ 2;
t565 = cos(qJ(2,1));
t625 = (pkin(1) * t565 + pkin(2) * t564) / t564 ^ 2;
t519 = sin(t531);
t542 = 0.1e1 / t555;
t624 = t519 * t542;
t520 = sin(t532);
t543 = 0.1e1 / t557;
t623 = t520 * t543;
t521 = sin(t533);
t544 = 0.1e1 / t559;
t622 = t521 * t544;
t522 = cos(t531);
t621 = t522 * t542;
t523 = cos(t532);
t620 = t523 * t543;
t524 = cos(t533);
t619 = t524 * t544;
t618 = t542 * t554;
t617 = t543 * t556;
t616 = t544 * t558;
t615 = t545 * t577;
t614 = t547 * t577;
t613 = t549 * t577;
t612 = t575 + t576;
t610 = t489 * t630;
t609 = t490 * t629;
t608 = t491 * t628;
t607 = t492 * t630;
t606 = t493 * t629;
t605 = t494 * t628;
t535 = -rSges(3,2) * t641 + Icges(3,6);
t495 = t535 * t560 - t536 * t554;
t604 = t495 * t508 * t545;
t496 = t535 * t562 - t536 * t556;
t603 = t496 * t509 * t547;
t497 = t535 * t564 - t536 * t558;
t602 = t497 * t510 * t549;
t601 = t554 * t627;
t600 = t577 * t627;
t599 = t556 * t626;
t598 = t577 * t626;
t597 = t558 * t625;
t596 = t577 * t625;
t595 = t545 * t618;
t579 = 1 / pkin(1);
t594 = t579 * t618;
t593 = t547 * t617;
t592 = t579 * t617;
t591 = t549 * t616;
t590 = t579 * t616;
t589 = t495 * t615;
t588 = t496 * t614;
t587 = t497 * t613;
t586 = 0.2e1 * rSges(3,3) ^ 2 + t612;
t585 = Icges(2,3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + t651 + Icges(3,1) / 0.2e1;
t584 = t586 * t643 + t585;
t534 = m(2) * rSges(2,2) - t641;
t569 = m(2) * rSges(2,1);
t583 = (-(-t569 + (-rSges(3,1) * t560 + rSges(3,2) * t554) * m(3)) * t561 - t534 * t555) * pkin(1);
t582 = (-(-t569 + (-rSges(3,1) * t562 + rSges(3,2) * t556) * m(3)) * t563 - t534 * t557) * pkin(1);
t581 = (-(-t569 + (-rSges(3,1) * t564 + rSges(3,2) * t558) * m(3)) * t565 - t534 * t559) * pkin(1);
t483 = t584 + t650;
t484 = t584 + t649;
t485 = t584 + t648;
t578 = pkin(1) ^ 2;
t580 = Icges(1,3) + ((2 * t578) + t586) * t643 + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + m(2) * t578 + t585;
t482 = t581 + t485;
t481 = t582 + t484;
t480 = t583 + t483;
t479 = t580 + 0.2e1 * t581 + t648;
t478 = t580 + 0.2e1 * t582 + t649;
t477 = t580 + 0.2e1 * t583 + t650;
t476 = t587 + (t482 * t549 - t485 * t596) * t590;
t475 = t588 + (t481 * t547 - t484 * t598) * t592;
t474 = t589 + (t480 * t545 - t483 * t600) * t594;
t473 = (t482 * t619 + t485 * t605) * t579;
t472 = (t481 * t620 + t484 * t606) * t579;
t471 = (t480 * t621 + t483 * t607) * t579;
t470 = (t482 * t622 + t485 * t608) * t579;
t469 = (t481 * t623 + t484 * t609) * t579;
t468 = (t480 * t624 + t483 * t610) * t579;
t467 = (t479 * t619 + t482 * t605) * t579;
t466 = (t478 * t620 + t481 * t606) * t579;
t465 = (t477 * t621 + t480 * t607) * t579;
t464 = (t479 * t622 + t482 * t608) * t579;
t463 = (t478 * t623 + t481 * t609) * t579;
t462 = (t477 * t624 + t480 * t610) * t579;
t461 = t488 * t613 + (t479 * t549 - t482 * t596) * t590;
t460 = t487 * t614 + (t478 * t547 - t481 * t598) * t592;
t459 = t486 * t615 + (t477 * t545 - t480 * t600) * t594;
t1 = [m(4) + (t465 * t621 + t466 * t620 + t467 * t619 + (t471 * t633 + t472 * t632 + t473 * t631) * t577) * t579, (t465 * t624 + t466 * t623 + t467 * t622 + (t471 * t636 + t472 * t635 + t473 * t634) * t577) * t579, (t465 * t595 + t466 * t593 + t467 * t591 + ((t492 * t604 + t493 * t603 + t494 * t602) * t577 + (-t473 * t597 + t524 * t637) * t544 + (-t472 * t599 + t523 * t638) * t543 + (-t471 * t601 + t522 * t639) * t542) * t577) * t579; (t462 * t621 + t463 * t620 + t464 * t619 + (t468 * t633 + t469 * t632 + t470 * t631) * t577) * t579, m(4) + (t462 * t624 + t463 * t623 + t464 * t622 + (t468 * t636 + t469 * t635 + t470 * t634) * t577) * t579, (t462 * t595 + t463 * t593 + t464 * t591 + ((t489 * t604 + t490 * t603 + t491 * t602) * t577 + (-t470 * t597 + t521 * t637) * t544 + (-t469 * t599 + t520 * t638) * t543 + (-t468 * t601 + t519 * t639) * t542) * t577) * t579; (t459 * t621 + t460 * t620 + t461 * t619 + (t474 * t633 + t475 * t632 + t476 * t631) * t577) * t579, (t459 * t624 + t460 * t623 + t461 * t622 + (t474 * t636 + t475 * t635 + t476 * t634) * t577) * t579, m(4) + (t459 * t595 + t460 * t593 + t461 * t591) * t579 + ((t644 + t645 + t646) * t577 * (t612 * m(3) + Icges(3,3)) + ((t488 * t644 + (-t476 - t587) * t625) * t616 + (t487 * t645 + (-t475 - t588) * t626) * t617 + (t486 * t646 + (-t474 - t589) * t627) * t618) * t579) * t577;];
MX  = t1;
