% Calculate Gravitation load for parallel robot
% P3PRRRR1G2A0
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
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
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
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:16
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRRR1G2A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G2A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G2A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G2A0_gravload_para_pf_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR1G2A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR1G2A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRRR1G2A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G2A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G2A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:16:25
% EndTime: 2020-03-09 21:16:25
% DurationCPUTime: 0.29s
% Computational Cost: add. (237->73), mult. (519->149), div. (54->10), fcn. (342->18), ass. (0->66)
t598 = sin(qJ(3,1));
t604 = cos(qJ(3,1));
t639 = -rSges(3,1) * t604 + rSges(3,2) * t598;
t596 = sin(qJ(3,2));
t602 = cos(qJ(3,2));
t638 = -rSges(3,1) * t602 + rSges(3,2) * t596;
t594 = sin(qJ(3,3));
t600 = cos(qJ(3,3));
t637 = -rSges(3,1) * t600 + rSges(3,2) * t594;
t636 = g(3) * rSges(3,3);
t591 = legFrame(3,2);
t575 = sin(t591);
t578 = cos(t591);
t572 = t578 * g(1) - t575 * g(2);
t585 = 0.1e1 / t600;
t569 = t575 * g(1) + t578 * g(2);
t595 = sin(qJ(2,3));
t601 = cos(qJ(2,3));
t614 = g(3) * t601 + t569 * t595;
t629 = ((rSges(3,1) * t572 + t614 * rSges(3,2)) * t600 + t594 * (t614 * rSges(3,1) - rSges(3,2) * t572)) * t585;
t592 = legFrame(2,2);
t576 = sin(t592);
t579 = cos(t592);
t573 = t579 * g(1) - t576 * g(2);
t587 = 0.1e1 / t602;
t570 = t576 * g(1) + t579 * g(2);
t597 = sin(qJ(2,2));
t603 = cos(qJ(2,2));
t613 = g(3) * t603 + t570 * t597;
t628 = ((rSges(3,1) * t573 + t613 * rSges(3,2)) * t602 + t596 * (t613 * rSges(3,1) - rSges(3,2) * t573)) * t587;
t593 = legFrame(1,2);
t577 = sin(t593);
t580 = cos(t593);
t574 = t580 * g(1) - t577 * g(2);
t589 = 0.1e1 / t604;
t571 = t577 * g(1) + t580 * g(2);
t599 = sin(qJ(2,1));
t605 = cos(qJ(2,1));
t612 = g(3) * t605 + t571 * t599;
t627 = ((rSges(3,1) * t574 + t612 * rSges(3,2)) * t604 + t598 * (t612 * rSges(3,1) - rSges(3,2) * t574)) * t589;
t582 = 0.1e1 / t595;
t626 = t582 * t585;
t625 = t582 * t601;
t583 = 0.1e1 / t597;
t624 = t583 * t587;
t623 = t583 * t603;
t584 = 0.1e1 / t599;
t622 = t584 * t589;
t621 = t584 * t605;
t620 = t595 * t600;
t619 = t597 * t602;
t618 = t599 * t604;
t617 = t569 * t626;
t616 = t570 * t624;
t615 = t571 * t622;
t606 = rSges(2,2) * g(3);
t607 = rSges(2,1) * g(3);
t563 = ((-rSges(2,1) * t569 + t606) * t601 + (rSges(2,2) * t569 + t607) * t595) * m(2) + ((t637 * t569 - t636) * t601 + (-t569 * rSges(3,3) - t637 * g(3)) * t595) * m(3);
t611 = t563 / t600 ^ 2 * t594 * t625;
t564 = ((-rSges(2,1) * t570 + t606) * t603 + (rSges(2,2) * t570 + t607) * t597) * m(2) + ((t638 * t570 - t636) * t603 + (-t570 * rSges(3,3) - t638 * g(3)) * t597) * m(3);
t610 = t564 / t602 ^ 2 * t596 * t623;
t565 = ((-rSges(2,1) * t571 + t606) * t605 + (rSges(2,2) * t571 + t607) * t599) * m(2) + ((t639 * t571 - t636) * t605 + (-t571 * rSges(3,3) - t639 * g(3)) * t599) * m(3);
t609 = t565 / t604 ^ 2 * t598 * t621;
t608 = 0.1e1 / pkin(2);
t581 = m(1) + m(2) + m(3);
t1 = [-m(4) * g(1) + (-(t577 * t618 - t580 * t598) * t615 - (t576 * t619 - t579 * t596) * t616 - (t575 * t620 - t578 * t594) * t617) * t581 + (t578 * t611 + t579 * t610 + t580 * t609 + (-t578 * t629 - t579 * t628 - t580 * t627) * m(3)) * t608; -m(4) * g(2) + (-(t577 * t598 + t580 * t618) * t615 - (t576 * t596 + t579 * t619) * t616 - (t575 * t594 + t578 * t620) * t617) * t581 + (-t575 * t611 - t576 * t610 - t577 * t609 + (t575 * t629 + t576 * t628 + t577 * t627) * m(3)) * t608; -m(4) * g(3) + (-t563 * t626 - t564 * t624 - t565 * t622) * t608 + (-t569 * t625 - t570 * t623 - t571 * t621) * t581;];
taugX  = t1;
