% Calculate Gravitation load for parallel robot
% P3PRRR1G3P3A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2020-03-09 21:07
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRR1G3P3A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G3P3A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G3P3A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G3P3A0_gravload_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRR1G3P3A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRR1G3P3A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRR1G3P3A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G3P3A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G3P3A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:06:46
% EndTime: 2020-03-09 21:06:47
% DurationCPUTime: 0.25s
% Computational Cost: add. (483->73), mult. (480->118), div. (45->5), fcn. (309->18), ass. (0->63)
t589 = legFrame(3,2);
t578 = sin(t589);
t581 = cos(t589);
t560 = t581 * g(1) - t578 * g(2);
t584 = pkin(7) + qJ(2,3);
t575 = qJ(3,3) + t584;
t563 = sin(t575);
t566 = cos(t575);
t608 = rSges(3,2) * g(3);
t610 = rSges(3,1) * g(3);
t539 = (rSges(3,1) * t560 - t608) * t563 + t566 * (rSges(3,2) * t560 + t610);
t569 = sin(t584);
t572 = cos(t584);
t542 = 0.1e1 / (t563 * t572 - t569 * t566);
t616 = t539 * t542;
t590 = legFrame(2,2);
t579 = sin(t590);
t582 = cos(t590);
t561 = t582 * g(1) - t579 * g(2);
t585 = pkin(7) + qJ(2,2);
t576 = qJ(3,2) + t585;
t564 = sin(t576);
t567 = cos(t576);
t540 = (rSges(3,1) * t561 - t608) * t564 + t567 * (rSges(3,2) * t561 + t610);
t570 = sin(t585);
t573 = cos(t585);
t543 = 0.1e1 / (t564 * t573 - t570 * t567);
t615 = t540 * t543;
t591 = legFrame(1,2);
t580 = sin(t591);
t583 = cos(t591);
t562 = t583 * g(1) - t580 * g(2);
t586 = pkin(7) + qJ(2,1);
t577 = qJ(3,1) + t586;
t565 = sin(t577);
t568 = cos(t577);
t541 = (rSges(3,1) * t562 - t608) * t565 + t568 * (rSges(3,2) * t562 + t610);
t571 = sin(t586);
t574 = cos(t586);
t544 = 0.1e1 / (t565 * t574 - t571 * t568);
t614 = t541 * t544;
t613 = (pkin(2) * t572 + pkin(3) * t566) * t616;
t612 = (pkin(2) * t573 + pkin(3) * t567) * t615;
t611 = (pkin(2) * t574 + pkin(3) * t568) * t614;
t609 = rSges(2,2) * g(3);
t607 = m(3) / pkin(3);
t606 = pkin(2) * t569;
t605 = pkin(2) * t570;
t604 = pkin(2) * t571;
t588 = g(3) * pkin(2) * m(3);
t592 = rSges(2,1) * g(3);
t600 = t542 * ((t588 + m(2) * (rSges(2,2) * t560 + t592)) * t572 + m(2) * (rSges(2,1) * t560 - t609) * t569 + (t560 * t606 + t539) * m(3));
t599 = t543 * ((t588 + m(2) * (rSges(2,2) * t561 + t592)) * t573 + m(2) * (rSges(2,1) * t561 - t609) * t570 + (t561 * t605 + t540) * m(3));
t598 = t544 * ((t588 + m(2) * (rSges(2,2) * t562 + t592)) * t574 + m(2) * (rSges(2,1) * t562 - t609) * t571 + (t562 * t604 + t541) * m(3));
t597 = t566 * t600;
t596 = t567 * t599;
t595 = t568 * t598;
t594 = 0.1e1 / pkin(2);
t587 = m(1) + m(2) + m(3);
t559 = t580 * g(1) + t583 * g(2);
t558 = t579 * g(1) + t582 * g(2);
t557 = t578 * g(1) + t581 * g(2);
t1 = [-m(4) * g(1) + (-t578 * t557 - t579 * t558 - t580 * t559) * t587 + (t581 * t597 + t582 * t596 + t583 * t595 + (-t581 * t613 - t582 * t612 - t583 * t611) * t607) * t594; -m(4) * g(2) + (-t557 * t581 - t558 * t582 - t559 * t583) * t587 + (-t578 * t597 - t579 * t596 - t580 * t595 + (t578 * t613 + t579 * t612 + t580 * t611) * t607) * t594; -m(4) * g(3) + (-t563 * t600 - t564 * t599 - t565 * t598 + ((pkin(3) * t565 + t604) * t614 + (pkin(3) * t564 + t605) * t615 + (pkin(3) * t563 + t606) * t616) * t607) * t594;];
taugX  = t1;
